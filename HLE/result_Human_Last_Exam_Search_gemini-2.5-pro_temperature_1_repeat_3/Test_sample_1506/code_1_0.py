def primes(n, is_prime_flag):
  prime_list = []
  for i in range(n + 1):
    if is_prime_flag[i]:
      prime_list.append(i)
  return prime_list