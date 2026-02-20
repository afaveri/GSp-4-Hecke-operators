#!/usr/bin/env python
# coding: utf-8

# In[2]:


var('p')

# Here q = p^(1/2) 
LR_q.<q> = LaurentPolynomialRing(QQ, 1)
LR_q_yz.<y,z> = LaurentPolynomialRing(LR_q, 2)

# Ambient ring (to avoid nested polynomial rings in computations)
LR.<q,y,z> = LaurentPolynomialRing(QQ, 3)


# In[ ]:





# In[ ]:


############################## Computation of volumes ##############################


# In[ ]:





# In[3]:


# In what follows, t = t(n, m, l) is the Hecke operator corresponding to diag(p^n, p^m, p^(l-n), p^(l-m))
# Choose positive roots a_1(t) = m - n, a_2(t) = l - 2*n, a_3(t) = l - n - m, a_4(t) = l - 2*m
# Thus the half-sum of positive roots is rho(t) = (3/2)*l - m - 2*n
# With this choice, t is positive if and only if n <= m <= l/2

# Weyl group elements

def w1(n,m,l): # (1)(2)(3)(4)
    return (n,m,l)

def w2(n,m,l): # (12)(3)(4)
    return (m,n,l)

def w3(n,m,l): # (13)(2)(4)
    return (l-n,m,l)

def w4(n,m,l): # (1)(24)(3)
    return (n,l-m,l)

def w5(n,m,l): # (13)(24)
    return (l-n,l-m,l)

def w6(n,m,l): # (1432)
    return (m,l-n,l)

def w7(n,m,l): # (1234)
    return (l-m,n,l)

def w8(n,m,l): # (14)(23)
    return (l-m,l-n,l)

Weyl_group = {
    1: w1,
    2: w2,
    3: w3,
    4: w4,
    5: w5,
    6: w6,
    7: w7,
    8: w8,
}

def w(j, n, m, l):
    try:
        return Weyl_group[j](n, m, l)
    except KeyError:
        raise ValueError(f"Invalid Weyl group index j={j}. Must be between 1 and 8.")

# Compute <t, rho>.
def rho_inner(n,m,l):
    return (3/2)*l - m - 2*n

# Compute sup_{w in W} <w.t, rho>.
def G_norm(n,m,l):
    res = 0
    for j in range(1,9):
        (a,b,c) = w(j, n, m, l)
        res = max(res, rho_inner(a,b,c))  
    return res

# Let C_0 = {t(n,m,l) : n < m < l/2}, and for w in W let C_0(w) = {1 <= i <= 4 : a_i(w.C) < 0 for all C in C_0}.
C0_data = {
    1: [],
    2: [1],
    3: [1, 2, 3],
    4: [4],
    5: [1, 2, 3, 4],
    6: [3, 4],
    7: [1, 2],
    8: [2, 3, 4],
}

# Return a list of 1 <= i <= 4 such that i is in C_0(w).
def C_0(j):
    try:
        return C0_data[j]
    except KeyError:
        raise ValueError(f"Invalid Weyl group index j={j}. Must be between 1 and 8.")

# Return W_t = {w in W : w.t = t}, for t = t(n,m,l).
def W_t(n,m,l):
    return [j for j in range(1, 9) if (n, m, l) == w(j, n, m, l)]

# Compute Q_t(1/p) = sum_{w in W_t} p^(- |C_0(w)|) for t = t(n,m,l).
def Q_t(n,m,l):
    return sum(q^(- 2 * len(C_0(j))) for j in W_t(n,m,l))

# Compute Q(1/p) = sum_{w in W} p^(- |C_0(w)|).
def Q():
    return sum(q^(- 2 * len(C_0(j))) for j in range(1, 9))

# Compute vol(K_p t K_p), with vol(K_p) = 1. That is, compute the number of single cosets K_p t K_p / K_p.
# It is equal to p^(2 * || t ||^*) * Q(1/p) / Q_t(1/p). This only works for positive t.
def q_volume(n,m,l):
    Q_factor = Q() / Q_t(n,m,l)
    return Q_factor * q^(4 * G_norm(n,m,l))

# This only works for positive t.
def volume(n,m,l):
    return q_volume(n,m,l).subs(q=sqrt(p))


# In[ ]:





# In[12]:


############################## Computation of roots ##############################


# In[ ]:





# In[13]:


# Positive roots

def a1(n,m,l):
    return m - n

def a2(n,m,l):
    return l - 2*n

def a3(n,m,l):
    return l - n - m

def a4(n,m,l):
    return l - 2*m

positive_roots = {
    1: a1,
    2: a2,
    3: a3,
    4: a4,
}

def pos_root(i, n, m, l):
    try:
        return positive_roots[i](n, m, l)
    except KeyError:
        raise ValueError(f"Invalid positive root index i={i}. Must be between 1 and 4.")

# Returns True if and only if n <= m <= l/2.
def is_positive(n,m,l):
    for i in range(1,5):
        if pos_root(i, n, m, l) < 0:
            return False
    return True


# In[ ]:





# In[19]:


############################## Macdonald's formula ##############################


# In[ ]:





# In[20]:


# Denote t = t(n, m, l)

# Computes s(t).
def s_poly(n,m,l):
    return y^(m-n) * z^(l-2*n)

# Computes s(w.t) = (w^{-1}.s)(t).
def s_poly_Weyl(j, n, m, l):
    (d,e,f) = w(j, n, m, l)
    return s_poly(d,e,f)

# Computes (1 - s(w.t)^{-1} p^{-1}) / (1 - s(w.t)^{-1}).
def c_factor(j, n, m, l):
    return (1 - s_poly_Weyl(j, -n, -m, -l) / q^2) / (1 - s_poly_Weyl(j, -n, -m, -l))

# Computes c(w^{-1}.s).
def c_final(j):
    return c_factor(j, -1, 1, 0) * c_factor(j, -1, 0, 0) * c_factor(j, -1, -1, 0) * c_factor(j, 0, -1, 0)

# Computes sum_{w in W} (w.s)(t) * c(w.s).
def pol(n,m,l):
    return sum(s_poly_Weyl(j, n, m, l) * c_final(j) for j in range(1, 9))

# Spherical function omega_s(t). Assumes t is positive. Note that it contains factors Q(1/p)^(-1) and p^{ - || t ||^*}.
def spherical(n,m,l):
    return (Q())^(-1) * q^(- 2 * G_norm(n,m,l)) * pol(n,m,l)

# Spherical transform of indicator of K_p t K_p. Assumes t is positive.
def satake(n,m,l):
    return q_volume(n,m,l) * spherical(n,m,l)


# In[ ]:





# In[21]:


############################## Decomposition into basic Hecke operators ##############################


# In[ ]:





# In[22]:


# Decompose Satake polynomial f into the linear basis { satake(0,m,l) } with 0 <= 2*m <= l <= lim.
# Returns a dictionary mapping (m,l) to its coefficient, or None if no solution is found (in that case one must increase lim).
def decompose_satake(f, lim=0):

    # Consider f as Laurent polynomial in y and z
    f = LR_q_yz(f)

    # Index set for basis elements
    index_set = [(m, l)
                 for l in range(lim + 1)
                 for m in range(floor(l/2) + 1)]

    basis = [LR_q_yz(satake(0, m, l)) for (m, l) in index_set]

    # Create list of all relevant monomials
    monomials = set(f.monomials())
    for b in basis:
        monomials |= set(b.monomials()) # set union

    monomials = list(monomials)

    # Build linear system: basis_coeffs * x = f_coeffs, over the Laurent polynomial ring LR_q
    R = LR_q # coefficient ring
    basis_coeffs = matrix(R, len(monomials), len(basis)) 
    f_coeffs = vector(R, len(monomials))

    for i, mon in enumerate(monomials):
        f_coeffs[i] = f.monomial_coefficient(mon)
        for j, basis_j in enumerate(basis):
            basis_coeffs[i, j] = basis_j.monomial_coefficient(mon)

    # Solve linear system
    try:
        sol = basis_coeffs.solve_right(f_coeffs)
    except: # Could fail if lim is not large enough
        return decompose_satake(f, lim + 1) # increase lim and try again

    # Verify exact reconstruction
    reconstruction = sum(LR_q_yz(coeff) * polyn for coeff, polyn in zip(sol, basis))
    if reconstruction != f: # Something is wrong
        print(f'Reconstruction failed!')
        return None

    # Return coefficients as a dictionary indexed by (m,l)
    Rat.<q> = PolynomialRing(QQ, 1) # cast coefficients to rational polynomial in q
    return {index : Rat(coeff) for coeff, index in zip(sol, index_set)}


# In[23]:


# Return absolute value of a univariate polynomial in q (for sufficiently large q)
def poly_abs(poly):
    if poly == 0:
        return 0
    elif poly.leading_coefficient() < 0:
        return -poly
    else:
        return poly

# Compute exact L1(G) norm from decomposition.
def compute_L1_G(solution):
    L1_G = 0

    for (m, l), coeff in solution.items():
        L1_G += poly_abs(coeff) * q_volume(0, m, l)

    return L1_G

# Compute the approximate L1(H) norm, keeping only leading power in p (correct up to constant multiple).
def estimate_L1_H(solution):
    L1_H = 0

    for (m, l), coeff in solution.items():
        if m == 0:
            L1_H += poly_abs(coeff) * q^(4*l) # approximate contribution (up to constants) to L1(H) -- comes only from satake(0, 0, l)

    if L1_H == 0:
        return 0

    # Since all bounds are up to upper and lower multiplicative constants, can return only leading term
    return q^(L1_H.degree())


# In[24]:


# Display decomposition of Satake polynomial f into basic Hecke operators, where (m, l) corresponds to coefficient of t(0, m, l).
def find_sols(f, lim=0):

    solution = decompose_satake(f, lim)

    if solution is None:
        print(f'No solution found, something went wrong!')
        return

    L1_G = compute_L1_G(solution)
    L1_H = estimate_L1_H(solution)

    p_solution = {index: coeff.subs(q=sqrt(p)) for index, coeff in solution.items()}
    L1_G = L1_G.subs(q=sqrt(p))
    L1_H = L1_H.subs(q=sqrt(p))

    print(f'L1(G) norm = {L1_G}')
    print(f'L1(H) norm << {L1_H}')

    return p_solution


# In[ ]:





# In[25]:


############################## Some tests ##############################


# In[ ]:





# In[26]:


T1 = satake(0, 0, 1)
T2 = satake(0, 1, 2)
T3 = satake(0, 0, 2)


# In[27]:


find_sols(1)


# In[28]:


find_sols(T1)


# In[29]:


find_sols(T2)


# In[30]:


find_sols(T3)


# In[ ]:





# In[31]:


find_sols(T1 - T2) # check that L1(G) norm works (should be around p^4, not 0)


# In[32]:


find_sols(T1^2) # Check relation T1^2 − (p + 1)*T2 − T3 = 1 + p + p^2 + p^3, which is (17) of our non-escape of mass paper


# In[33]:


# Check invariance under center
(n,m,l) = (0, 2, 7)
cent = 5
(d,e,f) = center(cent, n, m, l)
satake(n,m,l) == satake(d,e,f)


# In[ ]:





# In[34]:


############################## Operators used in Lemma 9 and their decompositions ##############################


# In[ ]:





# In[35]:


find_sols(T2)


# In[36]:


find_sols(T2^2)


# In[ ]:





# In[37]:


sigma_p = T2^2 - (q^2+1) * T1^2 # Second operator from paper 


# In[38]:


find_sols(sigma_p) # Get cancellation on L1(H) norm (from p^5 to p^3)


# In[39]:


find_sols(sigma_p^2) # Get L1(H) norm p^10, so no cancellation (but will gain on these diagonal terms from amplification)


# In[ ]:




