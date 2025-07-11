def solve_vtable_loads():
    """
    Analyzes a C++ snippet to determine the number of vtable loads
    under perfect compiler optimization.
    """

    # Call 1: a->foo()
    # The compiler knows the exact type of 'a' is 'A' from 'new A()'.
    # It can devirtualize the call, making it a direct call to A::foo().
    # No vtable lookup is needed.
    call1_loads = 0

    # Call 2: a->foo() after escape(a)
    # The 'escape(a)' function makes the dynamic type of 'a' unknown
    # to the compiler. A full virtual dispatch is required.
    # This requires loading the vptr to find the vtable.
    call2_loads = 1

    # Call 3: b->foo() after new(a) B
    # The placement new 'new(a) B' tells the compiler that the
    # object's dynamic type is now exactly 'B'.
    # The call can be devirtualized to a direct call to B::foo().
    # No vtable lookup is needed.
    call3_loads = 0

    total_loads = call1_loads + call2_loads + call3_loads

    print("Analysis of virtual table loads with perfect optimizations:")
    print("-" * 50)
    print(f"First call 'a->foo()': The compiler knows the object's type is 'A', so it can devirtualize the call. Loads: {call1_loads}")
    print(f"Second call 'a->foo()': After 'escape(a)', the type is unknown, requiring a true virtual dispatch. Loads: {call2_loads}")
    print(f"Third call 'b->foo()': After placement 'new ... B', the compiler knows the type is 'B' and devirtualizes. Loads: {call3_loads}")
    print("-" * 50)
    print("Total virtual table loads:")
    # Final equation as requested
    print(f"{call1_loads} + {call2_loads} + {call3_loads} = {total_loads}")

solve_vtable_loads()