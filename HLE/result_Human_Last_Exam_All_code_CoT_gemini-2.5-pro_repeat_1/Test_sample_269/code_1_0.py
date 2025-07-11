def simulate_cosmo_loco():
    """
    This script simulates the execution of the provided Scheme code,
    which uses `call/cc` (call-with-current-continuation) to create
    non-linear control flow.

    The simulation uses a state machine approach to track the program's
    location and the state of the continuation variable 'k'.
    """
    
    # --- State Variables ---
    # pc: The "Program Counter", indicating where we are in the code.
    # k:  Stores the current continuation function.
    # cc_return_val: Holds the value passed when a continuation is invoked.
    # output: A list to store the numbers printed.
    # limit: A safety break to prevent a true infinite loop in the simulation.
    
    pc = 'start_let'
    k = None
    cc_return_val = None
    output = []
    limit = 12

    while len(output) < limit:
        
        # State 1: Start of the `let` block. Handles the initial `call/cc`.
        if pc == 'start_let':
            # This function represents the continuation for the `let` block itself.
            # When called, it jumps execution back here.
            def cont0(val):
                nonlocal pc, cc_return_val
                pc = 'start_let'
                cc_return_val = val

            if cc_return_val is None:
                # First time through: (call/cc) returns the new continuation.
                k = cont0
            else:
                # Arrived via a jump: (call/cc) returns the value it was called with.
                k = cc_return_val

            cc_return_val = None # Reset for the next call/cc
            pc = 'after_let' # Proceed to the `let` body

        # State 2: After binding k, execute `(write 1)`.
        elif pc == 'after_let':
            output.append(1)
            pc = 'after_w1_callcc'

        # State 3: The `call/cc` after `(write 1)`.
        elif pc == 'after_w1_callcc':
            # This continuation represents the code starting from `(write 2)`.
            def cont_after_w1(unused_val):
                nonlocal pc
                pc = 'after_w2'
            
            # Execute `(k c)`, where c is `cont_after_w1`.
            k(cont_after_w1)

        # State 4: Execute `(write 2)`.
        elif pc == 'after_w2':
            output.append(2)
            pc = 'after_w2_callcc'

        # State 5: The `call/cc` after `(write 2)`. This is where the loop happens.
        elif pc == 'after_w2_callcc':
            # This continuation represents the code starting from `(write 3)`.
            def cont_after_w2(unused_val):
                nonlocal pc
                pc = 'after_w3'
            
            # Execute `(k c)`. 'k' is still the continuation that points
            # to `after_w2`, so this jumps back, creating the infinite loop.
            k(cont_after_w2)

        # The following states are defined but will not be reached due to the loop.
        elif pc == 'after_w3':
            output.append(3)
            pc = 'after_w3_callcc'
            
        elif pc == 'after_w3_callcc':
            def cont_after_w3(unused_val):
                nonlocal pc
                pc = 'after_w4'
            k(cont_after_w3)
            
        elif pc == 'after_w4':
            output.append(4)
            pc = 'end'
        
        else: # pc == 'end'
            break
            
    # Print the output captured during the simulation.
    print(' '.join(map(str, output)), '...')

simulate_cosmo_loco()