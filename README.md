# GSp(4) Hecke operators

This is a companion to the paper

> ***Mass equidistribution for lifts on hyperbolic 4-manifolds***, by *Alexandre de Faveri* and *Zvi Shem-Tov*.

It provides a Sage implementation of the Satake isomorphism for the Hecke algebra of $GSp_4(\mathbb{Q}\_{p})$.

## Contents

The Jupyter notebook `GSp(4).ipynb` provides two main functions: `satake(n, m, l)`, which computes the image under the Satake isomorphism of (the double coset corresponding to) $T_{n,m,l} := \mathrm{diag}(p^n,p^m,p^{l-n},p^{l-m})$, and `find_sols(f)`, which decomposes the inverse image of an arbitrary $f$ into a linear combination of double cosets.

Actually we work modulo the center, i.e. over $PGSp_4(\mathbb{Q}\_{p})$, so this decomposition is in terms of (the double cosets corresponding to) $T\_{0,m,l}$ satisfying the positivity condition $0 \leq 2m \leq l$. The code can be easily adapted to handle the general case.

### Mathematical conventions

See Appendix B of the aforementioned paper for the conventions used in the code.

## Running

Start Sage’s Jupyter notebook

```bash
sage -n jupyter
```

and open `GSp(4).ipynb`.

### Requirements

- Tested with SageMath 10.8

## Acknowledgements

Thanks to Elad Zelingher for his help with SageMath.
