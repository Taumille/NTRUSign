## NTRUSign 

This repository creates an implementation of NTRUSign algorithm based on the paper [_NTRUSign: Digital Signatures Using the NTRU Lattice_](https://www.math.brown.edu/jpipher/NTRUSign_RSA.pdf) from Jeffrey Hopstein, Nick Howgrave-Graham, Jill Pipher, Joseph H. Silverman and William Whyte.
This implementation was made in Python with the packages :
- Numpy
- hashlib
- random

### Architecture

**Polynome.py** : This file contains every functions needed to deal with polynomials. It contains the class polynomials with the classical operations ($\times$, $+$, $-$, $\div$, modulo) but also specials functions, like polynomial inversion modulo $q$ and the "_star multiplication_".

#### Polynome.py

##### Star multiplication
The star multiplication, denoted by the symbol $*$ is an operation on the ring $R=\mathbb{Z}[X]/(X^N+1)$ which can be defined like this :

Let $(p,q) \in R^2$, with respectively $p_i$ and $q_i$ there i-th coefficient, if $r=p*q$ then
$\displaystyle r_k = \sum_{i+j\equiv k [N]} f_i\cdot g_j$

As $p$ and $q$ are of degree N, we can rewrite this equation: 
$$\displaystyle r_k = \sum_{i = 0}^{k}f_i\cdot g_{k-i} + \sum_{i = k+1}^{N}f_i\cdot g_{N+k-i}$$

Which can be easily computed

##### Inverse mod q 

This method find the inverse of a polynomial taking the star multiplication as the multiplication method, so let $p\in R$ and $q$ power of a prime, if $r$ is the inverse of p mod q, then 
$$p*r\equiv 1 [q]$$

N.B. : The modular operation in $R$ is calculated by applying this modulant to every coefficient.


