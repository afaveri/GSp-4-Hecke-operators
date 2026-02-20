# GSp(4) Hecke operators

We give a Sage implementation of the Satake isomorphism for the group $GSp_4(\mathbb{Q}_p)$.

The Jupyter notebook GSp(4).ipynb provides two main functions: satake(n, m, l), which computes the image of (the double coset corresponding to) 
$T_{n, m, l} := \mathrm{diag}(p^n, p^m, p^{l-n}, p^{l-m})$, and find_sols(f), which decomposes the pre-image of $f$ into double cosets. 

Actually we work modulo the center, i.e. in $PGSp_4(\mathbb{Q}_p)$, so this decomposition is in terms of (the double cosets corresponding to) $T_{0, m, l}$ 
satisfying the positivity condition $0 \leq m \leq l/2$. The code can be easily adapted to handle the general case.

## Acknowledgements

Thanks to Elad Zelingher for his help with SageMath.
