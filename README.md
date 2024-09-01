# thesis-tempname
An ISD improvement for LDPC codes

The repo follows the following structure:

- H-generation: contains scripts for structured parity check matrix generation
- costs: contains scripts for computational complexity estimation of involved algorithms
- isds: contians different isd implementations in both python (and C in the near future) 
- tex: contains the actual graduation thesis material (trying to understand how to link it via overleaf)
- references: contains some downloaded papers I have read regarding ISD and codes

## H-generation

## Isds

Some comments on the files contained in isds directory:

- isd: contains first stern and lee and brickell implementation (finds all codewords)
- pge: code for PGE execution
- sparsestern: contains first sparsestern implementation (finds all codewords)
- decisional\_sparsestern1: sparsestern implementation that stops when a codeword is found (mutated from sparsestern)
- decisional\_sparsestern2: sparsestern implementation that doesn't use total list enumeration for MITM (mutated from above version)
- decisional\_sparsestern3: sparsestern implementation that doesn't shift A matrix right (mutated from above version)
- milp\_sparsestern: sparsestern implementation that uses a solver (SCIP) to get optimal z (mutuated from sparsestern)
- decisional\_stern: newest stern implementation (mutated from decisional\_sparsestern3)

## Costs

Some comments on the files contained in isds directory:

- costs_idea: comparison for fixed column weight
- isd_improvement:
- min_distance:
- Stern_small_subcodes_setting:
- Stern_small_subcodes_setting_approx:

## Tex

Contains current thesis progress, thus it's mainly latex files. Only updated from time to time.
I'm planning on building my own class. Still making my mind up about it.

## References

Currently contained references:

-


