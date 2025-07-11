def simulate_make():
    """
    This function simulates the `make all` command based on the provided Makefile
    and initial file state. It prints the final list of files in the directory.
    """

    # 1. Initial files and their timestamps.
    # A higher number means a newer file.
    files = {
        'X': 1051,
        'Y': 1052,
        'Z': 1054,
        'OPPS': 1132,
        'Makefile': 1134,
    }

    # Use a counter for timestamps of newly created files
    time_counter = 2000

    # Memoization to track processed targets to avoid re-evaluation
    processed_targets = set()
    
    # Track the call stack to detect circular dependencies
    call_stack = []

    def get_timestamp(filename):
        return files.get(filename, -1)

    def touch(filename):
        nonlocal time_counter
        files[filename] = time_counter
        time_counter += 1

    def build_target(target):
        if target in processed_targets:
            return

        if target in call_stack:
            # Handle circular dependency
            return
            
        call_stack.append(target)

        # Define Makefile rules inside the function scope
        rules = {
            'all': {'deps': ['T', 'Z', 'X', 'Opps'], 'cmd': None},
            'T': {'deps': ['Opps', 'X'], 'cmd': lambda: touch('A')},
            'X': {'deps': ['Y'], 'cmd': lambda: touch('B')},
            'Z': {'deps': ['Y'], 'cmd': lambda: touch('C')},
            'Opps': {'deps': ['T', 'Z'], 'cmd': lambda: touch('T')}
        }

        # If a target has no rule, it's a source file.
        if target not in rules:
            processed_targets.add(target)
            call_stack.pop()
            return
        
        rule = rules[target]
        dependencies = list(rule['deps']) # Make a mutable copy

        # Special handling for the known circular dependency `Opps <- T`
        if target == 'Opps' and 'T' in call_stack:
            dependencies.remove('T')

        # Recursively build dependencies
        for dep in dependencies:
            build_target(dep)

        # Check if the target needs to be rebuilt
        needs_build = False
        target_timestamp = get_timestamp(target)

        if target_timestamp == -1: # Target file doesn't exist
            needs_build = True
        else:
            for dep in dependencies:
                dep_timestamp = get_timestamp(dep)
                # Rebuild if any dependency is newer
                if dep_timestamp > target_timestamp:
                    needs_build = True
                    break
        
        # Execute command if needed
        if needs_build and rule['cmd']:
            rule['cmd']()
        
        processed_targets.add(target)
        call_stack.pop()

    # Start the simulation with the 'all' target
    build_target('all')
    
    # Get the final list of files and print them
    final_files = sorted(files.keys())
    print(", ".join(final_files))

# Run the simulation
simulate_make()