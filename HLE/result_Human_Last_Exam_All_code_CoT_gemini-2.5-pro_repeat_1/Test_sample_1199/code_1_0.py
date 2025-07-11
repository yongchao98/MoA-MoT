# Let's break down each virtual function call.
# A "vtable load" means the program has to read the object's vtable pointer
# to find the virtual function table and resolve the call at runtime.
# Compilers can optimize this away if the object's exact type is known at compile time.
# This optimization is called "devirtualization".

# Call 1: The first a->foo()
# Right after 'A* a = new A()', the compiler knows for certain that 'a' points to an
# object of type 'A'. With perfect optimization, it will devirtualize the call.
# It becomes a direct call to A::foo(), not a virtual one.
call_1_loads = 0

# Call 2: The second a->foo()
# This call happens after 'escape(a)'. The 'escape' function acts as an
# optimization barrier. The compiler can no longer know the real type of the
# object 'a' points to. It must perform a true virtual dispatch.
# This requires loading the vtable pointer from the object 'a'.
call_2_loads = 1

# Call 3: The b->foo()
# This call happens right after 'A* b = new(a) B'. The compiler knows from this
# statement that 'b' now points to an object of the exact type 'B'.
# So, it can devirtualize this call into a direct call to B::foo().
call_3_loads = 0

# Total vtable loads
total_loads = call_1_loads + call_2_loads + call_3_loads

print(f"Analysis of vtable loads with perfect compiler optimization:")
print(f"1. First call to a->foo(): Devirtualized, loads = {call_1_loads}")
print(f"2. Second call to a->foo() after escape(): True virtual dispatch, loads = {call_2_loads}")
print(f"3. Call to b->foo() after placement new: Devirtualized, loads = {call_3_loads}")
print(f"Total vtable loads = {call_1_loads} + {call_2_loads} + {call_3_loads} = {total_loads}")