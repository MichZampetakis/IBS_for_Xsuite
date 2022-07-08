import sys
import json
import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp
import matplotlib.pyplot as plt

from cpymad.madx import Madx
from scipy.constants import e as qe
from scipy.constants import m_p, m_e
#from lib.Xsuite_eval_emit_sig import *
#from lib.general_functions import *
from lib.IBSfunctions import *

## Read a MAD-X Sequence
n_slice_per_element = 4
mad = Madx()
mad.call('./chrom-corr_DR.newlattice_2GHz.seq')
mad.input(f'''
beam, particle=positron, energy=2.86, bunched;
use, sequence=RING;''')
twthick = mad.twiss().dframe()
mad.input(f'''
select, flag=MAKETHIN, SLICE={n_slice_per_element}, thick=false;
select, flag=MAKETHIN, pattern=wig, slice=1;
MAKETHIN, SEQUENCE=ring, MAKEDIPEDGE=true;
use, sequence=RING;''')
twthin = mad.twiss().dframe()
line = xt.Line.from_madx_sequence(mad.sequence['RING'])

# Need to find if the aperture is changed in Xsuite!!!!!! #
# line.insert_element(idx = 0, # start of the ring
#                     element = xt.LimitRect(min_x = -0.04, max_x = 0.04, min_y = -0.04, max_y = 0.04),
#                     name = 'aperture')


## Choose a context
context = xo.ContextCpu()         # For CPU

## Transfer lattice on context and compile tracking code
tracker = xt.Tracker(_context=context, line=line)

cavs = [ee for ee in line.elements if isinstance(ee, xt.Cavity)]
for cc in cavs:
    cc.lag = 180
    #if cc.voltage > 0:
    #    cc.voltage = -cc.voltage

p0 = xp.Particles(mass0=xp.ELECTRON_MASS_EV, q0=1, p0c=2.86e9)
tw = tracker.twiss(particle_ref = p0)

bunch_intensity = 4.4e9
sigma_z = 1.58e-3
n_part = int(5e3)
nemitt_x = 5.6644e-07
nemitt_y = 3.7033e-09
particles = xp.generate_matched_gaussian_bunch(
         num_particles=n_part, total_intensity_particles=bunch_intensity,
         nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
         particle_ref=p0, tracker=tracker)

# ----- Initialize IBS -----
IBS = NagaitsevIBS()
IBS.set_beam_parameters(particles)
IBS.set_optic_functions(tw)

#emit = Evaluate_Sigma_Emit(10,20,10,20,20,20)
#emit.define_bins_width(particles, tw)

## Initialize dictionaries
n_turns = 1000
turn_by_turn = {}

#record_coord = ['x', 'px', 'y', 'py', 'zeta', 'delta', 'state', 'particle_id']
#for nn in record_coord:
#    turn_by_turn[nn] = np.zeros((n_turns, n_part), dtype = getattr(particles, nn).dtype)

record_emit = ['eps_x', 'eps_y', 'sig_delta', 'bl']
for nn in record_emit:
    turn_by_turn[nn] = np.zeros((n_turns), dtype = float)
print('hahaha')    
## Track (saving turn-by-turn data)
for i in range(n_turns):
    print('Turn = ', i)
    print('N_part = ',len(particles.x[particles.state > 0]))
    
#    for nn in record_coord:
#        turn_by_turn[nn][i,:] = getattr(particles, nn)

    sig_x = np.std(particles.x[particles.state > 0])
    sig_y = np.std(particles.y[particles.state > 0])
    sig_delta = np.std(particles.delta[particles.state > 0])
    turn_by_turn['bl'][i] = np.std(particles.zeta[particles.state > 0])
    turn_by_turn['sig_delta'][i] = sig_delta
    turn_by_turn['eps_x'][i] = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
    turn_by_turn['eps_y'][i] = sig_y**2 / tw['bety'][0]
    #turn_by_turn['eps_x'][i], turn_by_turn['eps_y'][i], turn_by_turn['sig_delta'][i], turn_by_turn['bl'][i] = emit.eval_emits(particles)
    
    if (i % 50 == 0):
        IBS.calculate_simple_kick(particles)
    IBS.apply_simple_kick(particles)

    # if (i % 50 == 0):
    #     IBS.calculate_kinetic_coefficients(particles)
    # IBS.apply_kinetic_kick(particles)
    
    tracker.track(particles)
#    rec_id = turn_by_turn['particle_id'].copy()
#    for nn in record_coord:
#        temp = turn_by_turn[nn].copy()
#        turn_by_turn[nn][:, rec_id] = temp

Emitt = []
Emitt.append(turn_by_turn['eps_x'])
Emitt.append(turn_by_turn['eps_y'])
Emitt.append(turn_by_turn['sig_delta'])
Emitt.append(turn_by_turn['bl'])
np.savetxt(r"./xsuite_simple.txt", np.array(Emitt).T.tolist(), fmt="%e")

sys.exit()

n_turns = 100
tracker.track(particles, num_turns=n_turns,
              turn_by_turn_monitor=True)

## Turn-by-turn data is available at:
tracker.record_last_track.x
tracker.record_last_track.px
sys.exit()

print(tracker.record_last_track.x)
print(tracker.line.element_names)
#print(tracker.line.elements[0])


## Track (saving turn-by-turn data)
# loop with 1 turn .. to evaluate IBS etc
# (later) or append 1 element in line to do the IBS kick

#rec_sig_x = np.zeros(n_turns, dtype = np.float)
#rec_sig_delta = np.zeros(n_turns, dtype = np.float)

'''
f, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,5))
ax1.plot(particles.x, particles.px, '.')
ax2.plot(particles.y, particles.py, '.')
ax3.plot(particles.zeta, particles.delta, '.')
plt.show()
'''

''' # Radiation if off in the multipoles
myline = line.to_dict()
for i in range(len(myline['element_names'])):
    if myline["elements"][i]["__class__"] == "Multipole":
        if myline["elements"][i]["radiation_flag"] != 0 : print('Radiation!')
'''

'''
# Particle at (0,0)x3 remains stable.
n_part = 1
particles = xp.Particles(_context=context, p0c=2.86e9,
                         x=0, px=0, y=0, py=0, zeta=0, delta=0, )
## Track (saving turn-by-turn data)
n_turns = 100
tracker.track(particles, num_turns=n_turns,
              turn_by_turn_monitor=True)

## Turn-by-turn data is available at:
tracker.record_last_track.x
tracker.record_last_track.px
'''
