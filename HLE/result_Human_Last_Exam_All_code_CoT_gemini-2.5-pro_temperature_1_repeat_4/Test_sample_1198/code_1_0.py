def solve():
    """
    This function provides the answers to the two questions asked by the user.
    
    Question 1: Is there any a>0 real number for that floor(a^n) = n mod 2 for every n>0 integer?
    The answer is Yes. Such a number can be constructed by defining a sequence of nested intervals for 'a'.
    At each step n, the condition imposes that a^n must lie in a certain interval. The intersection of these
    intervals can be shown to be non-empty.

    Question 2: Is there any a>0 real number for that floor(a^n) = n mod 3 for every n>0 integer?
    The answer is No. A constructive attempt, similar to the one for modulo 2, fails.
    By trying to define a sequence of intervals for 'a', we find that at some point the interval
    becomes empty. For instance, starting with floor(a) = 1, we are forced into a path of choices
    for floor(a^n) that leads to a contradiction at n=5, as no integer with the required residue
    mod 3 exists in the determined range for a^5. This indicates no such 'a' exists.
    """
    answer_mod_2 = "Yes"
    answer_mod_3 = "No"
    print(f"{answer_mod_2},{answer_mod_3}")

solve()