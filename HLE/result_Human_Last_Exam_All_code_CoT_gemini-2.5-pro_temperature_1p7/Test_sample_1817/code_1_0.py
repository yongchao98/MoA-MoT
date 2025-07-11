def solve():
    """
    Identifies the inappropriate implementations based on the Domain Model pattern.
    - A: Inappropriate (Fat Controller, Anemic Domain Model).
    - B: Inappropriate (God Object).
    - C: Appropriate (Rich Domain Model).
    - D: Inappropriate (Anemic Domain Model via Transaction Script).
    - E: Appropriate (Rich Domain Model with Domain Services).
    """
    inappropriate_options = ['A', 'B', 'D']
    inappropriate_options.sort()
    result = ",".join(inappropriate_options)
    print(result)

solve()