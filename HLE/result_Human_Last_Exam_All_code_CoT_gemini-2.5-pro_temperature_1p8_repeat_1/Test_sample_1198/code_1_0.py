def solve():
    """
    Solves the problem of existence of a real number 'a' for two different modular conditions.
    
    The problem asks:
    1. Is there any a>0 real number for that floor(a^n) % 2 == n % 2 for every n>0 integer?
    2. Is there any a>0 real number for that floor(a^n) % 3 == n % 3 for every n>0 integer?
    
    The thinking process outlined above leads to the following conclusions:
    1. For the modulo 2 case, such a number 'a' can be constructed. The answer is Yes.
    2. For the modulo 3 case, a contradiction arises, showing that no such 'a' exists. The answer is No.
    """
    answer_mod_2 = "Yes"
    answer_mod_3 = "No"
    
    print(f"{answer_mod_2},{answer_mod_3}")

solve()
<<<Yes,No>>>