# Example running with a protein and cofactors

This example should be run on a computer with a GPU and CUDA available for OpenMM simulation.
```
bash run-example.sh > run-example.log 2>&1
```

Note -- the calculation might take a very long time after this line if you don't have OpenEye installed:
```
INFO:openmmforcefields.generators.template_generators:Generating a residue template for [H][c]1[c]([H])[c]([H])[c]([C@]2([H])[N]([H])[c]3[c]([H])[c]([H])[c]([C]([F])([F])[F])[c]([H])[c]3[C@@]3([H])[O][C@@]([H])([C]([H])([H])[N]([H])[S](=[O])(=[O])[C]([H])([H])[C]([H])([H])[N+]([H])([H])[C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[C@]32[H])[c]([H])[c]1[H] using openff-2.0.0
```
This is expected -- it's just AmberTools generating AM1-BCC charges.