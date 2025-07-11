# Plan:
# 1. The first call `a->foo()` requires 1 vptr load and 1 vfunction load.
vptr_loads_1 = 1
vfunction_loads_1 = 1
print(f"Call a->foo(): {vptr_loads_1} vptr load, {vfunction_loads_1} vfunction load")

# 2. The escape(a) call invalidates any compiler assumptions about the object pointed to by 'a'.
#    The compiler cannot reuse the vptr.
#    So, the second call `a->bar()` requires its own loads.
vptr_loads_2 = 1
vfunction_loads_2 = 1
print(f"Call a->bar() after escape(): {vptr_loads_2} vptr load, {vfunction_loads_2} vfunction load")

# 3. std::launder(a) also acts as an optimization barrier, telling the compiler it must
#    re-evaluate the object at the given address.
#    So, the third call `b->foo()` also requires its own loads.
vptr_loads_3 = 1
vfunction_loads_3 = 1
print(f"Call b->foo() after launder(): {vptr_loads_3} vptr load, {vfunction_loads_3} vfunction load")

# 4. Calculate the total number of loads.
total_vptr_loads = vptr_loads_1 + vptr_loads_2 + vptr_loads_3
total_vfunction_loads = vfunction_loads_1 + vfunction_loads_2 + vfunction_loads_3

print(f"\nTotal vptr loads = {vptr_loads_1} + {vptr_loads_2} + {vptr_loads_3} = {total_vptr_loads}")
print(f"Total vfunction loads = {vfunction_loads_1} + {vfunction_loads_2} + {vfunction_loads_3} = {total_vfunction_loads}")
print(f"\nFinal Answer: {total_vptr_loads} vptr loads, {total_vfunction_loads} vfunction loads")