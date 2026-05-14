from amuse.community import *
from amuse.community.nbody7.interface import Nbody7
from amuse.ic.plummer import new_plummer_model
from pathlib import Path
import shutil
import numpy as np

NN = 150
converter = nbody_system.nbody_to_si(NN*0.7 | units.MSun, 1.0 | units.parsec)
cluster = new_plummer_model(NN, convert_nbody=converter)

#Check, display files in the current working directory
workdir = Path.cwd()
print("Current working directory:", str(workdir))

archive = workdir / "test"
archive.mkdir(exist_ok=True)
for pat in ("CMPACT", "NBSTAT", "OUT3*", "fort*", "STOP*", "run.out", "err.out"):
    for f in workdir.glob(pat):
        shutil.move(str(f), archive / f.name)
outfile = str(workdir / "run.out")
errfile = str(workdir / "err.out")
print("Output will be written to:", outfile)
print("Errors will be written to:", errfile)


inst = Nbody7(
    converter,
    mode='cpu',
    redirection="file",          # selects the file path of the dispatcher
    redirect_stdout_file=outfile,
    redirect_stderr_file=errfile,
)
#inst = Nbody7(converter, redirection='none')

inst.initialize_code()

print("ARCHAIN parameters: CLIGHT NBH IDIS")
print(inst.get_archain_params())
print("Setting new ARCHAIN parameters: ")
inst.set_archain_params(2.5e4, 0, 1)
print("New ARCHAIN parameters:")
print(inst.get_archain_params())

inst.parameters.ETAI = 0.01
inst.parameters.ETAR = 0.01
inst.parameters.ETAU = 0.1
inst.parameters.RMIN = 3.5e-04
inst.parameters.DTMIN = 3.5e-05
inst.parameters.NNBMAX = 32
inst.parameters.DTADJ = 0.5
inst.parameters.DELTAT = 0.5
inst.parameters.TCRIT = 1.0e+03
inst.parameters.QE = 1.0e-2
inst.set_kz(27,1)
inst.set_kz(28,3)
inst.commit_parameters()

inst.particles.add_particles(cluster)
inst.particles.name = np.arange(1, NN+1)
inst.particles.kstar = np.zeros(NN, dtype=int)
inst.commit_particles()
print(f"First particle NAME KSTAR (before evolution): {inst.particles[0].name} {inst.particles[0].kstar}")
print(f"Last particle NAME KSTAR (before evolution): {inst.particles[-1].name} {inst.particles[-1].kstar}")

# The first evolve_model() call lazy-invokes CALL NBODY6 which
# in turn runs SCALE / FPOLY0 / FPOLY2 and populates ZKIN, POT.
# Querying energies BEFORE this would return zeros (the COMMON
# block hasn't been touched yet) so we must establish the
# baseline AFTER an initial step.

t_unit = converter.to_si(1.0 | nbody_system.time)
inst.evolve_model(2.0 * t_unit)
print("1 Nb time =", inst.tstar.in_(units.Myr))
print("1 Nb vel  =", inst.vstar.in_(units.kms))

inst.evolve_model(2.5 * t_unit)
e0 = inst.kinetic_energy + inst.potential_energy
inst.evolve_model(3.0 * t_unit)
e1 = inst.kinetic_energy + inst.potential_energy

ttot = inst.get_time()
print(f"Time after evolution: {ttot/t_unit} Nb time or {ttot.in_(units.Myr)}")

print(f"First particle NAME KSTAR (after evolution): {inst.particles[0].name} {inst.particles[0].kstar}")
print(f"Last particle NAME KSTAR (after evolution): {inst.particles[-1].name} {inst.particles[-1].kstar}")

# Typically, relative enery error is expected to be
# under 1% for short evolutions on a smooth Plummer model.
rel_err = abs((e1 - e0) / e0)
print(f"Relative energy error: {rel_err:.2e}")
assert rel_err < 0.01, "Energy error exceeds 1%"

inst.cleanup_code()
inst.stop()
