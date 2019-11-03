## PHYS4061 LAB4

* Last: 3/11/2019
* optimisation with steepest descent (SD) & conjugate gradient (CG)
* Lennard-Jones potential (LJ): fcc Ar crys
* Tersoff potential (TE): dia Si crys

The folder contains the following:

* `lattice.c`: core codes (major stuff here)
* `tool.c`: small helper library
* `vector.c`: vector library
* `out/`: coord & potential log during optimisation (look here for results!)

Quick Lookup: search " [##] " for key changeable params in codes

## Crystals

* simple cubic (sc), body-centered (bc), face-centered (fcc), diamond (dia)
* here only func for fcc and dia are included
* resp: `Fcc()` and `Dia()`
* periodic boundary condition (pbc) assumed throughout, e.g. when calc nbor

## Definitions

* K = cutoff dist
* a = latt const
* others are documented at the beginning of `lattice.c`
* atoms params, e.g. coord, nbor dist, are all global

## Potentials

* units:
	- length: ang
	- energy: eV
* **Lennard-Jones**
	- nature: pair potential for sparse gas
	- test system: Ar fcc 555 ; K=8.5 ; a=5.27
	- ref: online src
* **Tersoff**
	- nature: three-body potential w/ local env dependence
	- test system: Si dia 333 ; K=999.0 ; a=5.43
	- ref: _New empirical approach for the structure and energy of covalent systems_. Tersoff, 1988.
* note: a such chosen that dV/da = 0, i.e. a minimises total potential
* chosen params are indicated in the respective potential functions
* potential functions calc potential for an indiv atom
* total potential is all indiv potentials summed
* indiv forces on an atom are calculated by perturbing and re-calc indiv potential change

## Optimisation

* **SD** and **CG**
* idea: "suitably" allowing atoms to descent (or slide) along forces
* here (generalised) force F = -∇V (steepest gradient)
* algo: positions R += _alpha_ * F
* _alpha_ depends on the update scheme
* convergence when F -> 0 and total potential no longer changes
* see `SD()` and `CG()` for details
* note: the optimisation func are general: may pass in an arbitrary potential func

## Outputs

* coord log: `log.xyz`
* potential log: screen

## Program

* simply do `./run.sh`
* outputs redirected to `run.log`

## Results

Whole process normally takes less than 30s to complete.

Forced to break if taking over 100 steps. But usually converge within 80 steps.

See coord/log files in out/.