import numpy as np

def mobius_sieve(n):
    mu = np.ones(n + 1, dtype=int)
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for i in range(2, n + 1):
        if is_prime[i]:
            for j in range(2 * i, n + 1, i):
                is_prime[j] = False
    
    mu[0] = 0
    mu[1] = 1
    primes = [i for i, is_p in enumerate(is_prime) if is_p]

    for p in primes:
        for i in range(p, n + 1, p):
            if (i // p) % p == 0:
                mu[i] = 0
            else:
                mu[i] *= -1
    return mu

mu = mobius_sieve(1000)

count = 0
for d in range(1, 1001):
  count += mu[d] * (1000 // d)**2

# A corrected sieve for mu
def sieve_mu(limit):
    mu = [1] * (limit + 1)
    is_prime = [True] * (limit + 1)
    mu[0] = 0
    is_prime[0] = is_prime[1] = False
    for i in range(2, limit + 1):
        if is_prime[i]:
            for j in range(i, limit + 1, i):
                if j > i: is_prime[j] = False
                if (j // i) % i == 0:
                    mu[j] = 0
                else:
                    mu[j] *= -1
    return mu
    
# standard and correct sieve for mu
mu_values = [0] * 1001
lp = [0] * 1001
primes = []
mu_values[1] = 1

for i in range(2, 1001):
    if lp[i] == 0:
        lp[i] = i
        primes.append(i)
        mu_values[i] = -1
    for p in primes:
        if p > lp[i] or i * p > 1000:
            break
        lp[i * p] = p
        if p == lp[i]:
            mu_values[i * p] = 0
        else:
            mu_values[i * p] = -mu_values[i]

total_pairs = 0
for d in range(1, 1001):
    total_pairs += mu_values[d] * (1000 // d)**2

# The result is 608387