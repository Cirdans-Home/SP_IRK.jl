# Running the experiments in the Toeplitz cluster

We run the experiments in the paper using the GPU Nodes of the Toeplitz cluster
located ad the Green Data Center of the University of Pisa. Experiments need
to be run under SLURM. The bash script in this folder are an example of how
this can be done, and can be readily adapted for other machines.

1. First you select the running option into the `exec.sh` script,
   namely what method you want to use, the number of time steps, etc.
2. Then you run `./genscript.sh exec.sh` which will create all the script for
   any feasible thread/stage combination
3. Finally you put everything into the execution queue by running `./execscript.sh`

If you want to modify everything to run on a different cluster/machine, you
just need to modify the configuration in the first lines of `genscript.sh`, i.e.,
```bash
stages=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
method="gauss"
order=1
n=100
```

> [!CAUTION]
> The `exec.sh` script has several hard-coded folder telling the system where the code
> is located on the Toeplitz machine, you have to adapt them to run your case.

> [!WARNING]
> In our experience to have reasonable performance you need to have a version of Julia
> which is compiled on the actual machine you intend to use it. The precompiled versions
> seem to undergo a significant slowdon.