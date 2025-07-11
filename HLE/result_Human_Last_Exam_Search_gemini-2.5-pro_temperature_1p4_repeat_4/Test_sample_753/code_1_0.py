from decimal import Decimal, getcontext

getcontext().prec = 100 # Set precision to 100 digits

n_factorial = Decimal(3628800)
p10 = Decimal(42)
c_lambda_size = Decimal(25200)

result = c_lambda_size * (n_factorial ** p10)