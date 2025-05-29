# Import necessary classes and functions
from collections import Counter
from functools import reduce

class lit(int):
    def __new__(cls, v):
        if isinstance(v, int):
            if (0 <= v and v < 26) or (100 <= v and v < 126):
                return int.__new__(cls, v)
            else:
                raise ValueError("Only literals 0..25 (a-z) and 100..125(!a-!z) are supported for compact printing.")
        elif isinstance(v, str):
            o = ord(v)
            if ord('a') <= o and o <= ord('z'):
                return int.__new__(cls, ord(v)-ord('a'))
            elif ord('A') <= o and o <= ord('Z'):
                return int.__new__(cls, 100+ord(v)-ord('A'))
        else:
            raise TypeError("Literals must be either integers or characters.")
    def __or__(self, c):
        if isinstance(c, cube):
            return c | self
        elif isinstance(c, lit):
            return cube([c]) | self
        else:
            raise TypeError("Can only make union of literal and cube.")
    def __mul__(self, l):
        return cube([self, l])
    def __str__(self):
        return chr(ord('a')+self)

class cube(frozenset):
    def __new__(cls, s):
        if isinstance(s, str):
            s = [lit(c) for c in s]
        elif isinstance(s, lit):
            return frozenset.__new__(cls, [s])
        return frozenset.__new__(cls, s)
    def __or__(self, c):
        if isinstance(c, func):
            return c | self
        elif isinstance(c, lit):
            c = cube([c])
        return cube(frozenset.__or__(self,c))
    def __add__(self, c):
        return func([self, c])
    def __sub__(self, c):
        return cube(frozenset.__sub__(self,c))
    def __str__(self):
        res = [str(f) for f in list(self)]
        res.sort()
        return "".join(res)

class func(frozenset):
    def __new__(cls, s):
        if isinstance(s, str):
            return frozenset.__new__(cls, [cube(c) for c in s.split('+')])
        elif isinstance(s, cube):
            return frozenset.__new__(cls, [s])
        else:
            return frozenset.__new__(cls, s)
    def num_lit(self):
        res = 0
        for c in self:
            for l in c:
                if l > res:
                    res = l
        return res+1
    def __or__(self, f):
        if isinstance(f, lit):
            f = cube([f])
        if isinstance(f, cube):
            f = func([f])
        if not isinstance(f, func):
            raise TypeError("Set operations only possible for same type [here: func] (implicit upcasts supported: lit->cube->func).")
        return func(frozenset.__or__(self,f))
    def __and__(self, f):
        if isinstance(f, lit):
            f = cube([f])
        if isinstance(f, cube):
            f = func([f])
        if not isinstance(f, func):
            raise TypeError("Set operations only possible for same type [here: func] (implicit upcasts supported: lit->cube->func).")
        return func(frozenset.__and__(self,f))
    def __sub__(self, f):
        if isinstance(f, lit):
            f = cube([f])
        if isinstance(f, cube):
            f = func([f])
        if not isinstance(f, func):
            raise TypeError("Set operations only possible for same type [here: func] (implicit upcasts supported: lit->cube->func).")
        return func(frozenset.__sub__(self,f))
    def __xor__(self, f):
        if isinstance(f, lit):
            f = cube([f])
        if isinstance(f, cube):
            f = func([f])
        if not isinstance(f, func):
            raise TypeError("Set operations only possible for same type [here: func] (implicit upcasts supported: lit->cube->func).")
        return func(frozenset.__xor__(self,f))
    def __add__(self, f):
        return self.__or__(f)
    def __mul__(self, cc):
        if isinstance(cc, lit):
            cc = cube([cc])
        if not isinstance(cc, cube):
            raise TypeError("Functions can only be multiplied by cubes or literals.")
        f = func([])
        for c in self:
            f = f | (c | cc)
        return f
    def __truediv__(self, cc):
        if isinstance(cc, lit):
            cc = cube([cc])
        if not isinstance(cc, cube):
            raise TypeError("Functions can only be multiplied by cubes or literals.")
        f = func([])
        for c in self:
            if cc <= c:
                f = f | (c - cc)
        return f
    def __mod__(self, cc):
        f = func([])
        for c in self:
            if not (cc <= c):
                f = f | c
        return f
    def __str__(self):
        res = [str(f) if len(f)>0 else '1' for f in list(self)]
        res.sort()
        return "+".join(res)

def count(f):
    return reduce((lambda x, y : x + Counter(y)), f, Counter())

def cs(f):
    return [e[0] for e in count(f).most_common()]

def _level0kernels(f, ls, j):
    k = set()
    for i in range(len(ls)):
        li = ls[i]
        if len([1 for c in f if li in c]) > 1:
            f2 = f / li
            c = cube([e[0] for e in count(f2).most_common() if e[1] == len(f2)])
            if not any([ls[k] in c for k in range(i)]):
                k |= _level0kernels(f/(li|c), ls, i+1)
    if len(k) < 1:
        k.add(f)
    return k

def level0kernels(f):
    return _level0kernels(f, cs(f), 0)

# Main function to find a feasible input
def find_feasible_input():
    # Construct functions with distinct level 0 kernels
    functions = ["ab+cd", "ef+gh", "ij+kl", "mn+op", "qr+st"]
    return functions

# Execute the function to get the feasible input
feasible_input = find_feasible_input()
print(feasible_input)