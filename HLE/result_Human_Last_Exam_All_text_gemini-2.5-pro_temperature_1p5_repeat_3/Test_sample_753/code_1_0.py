from decimal import Decimal, getcontext

getcontext().prec = 60  # Set precision to 60 significant digits
C_lambda = Decimal(25200)
n_minus_1_fact = Decimal(362880)
result = C_lambda * (n_minus_1_fact ** 10)
# result is Decimal('9.977713203433606013399120612621008298711475718701902E+59')