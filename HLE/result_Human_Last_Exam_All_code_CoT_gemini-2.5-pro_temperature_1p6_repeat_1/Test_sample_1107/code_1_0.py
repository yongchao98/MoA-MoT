import collections

def simulate_make_all():
    """
    This script simulates the 'make all' command based on the provided Makefile and file states.
    It determines which files will be present in the directory after the command completes.
    """
    
    # 1. Initial state of the directory with timestamps.
    # Note: The makefile refers to 'Opps', but the file is 'OPPS'.
    # On a case-sensitive system, these are different.
    files = {
        'X': 1051, 
        'Y': 1052, 
        'Z': 1054, 
        'OPPS': 1132, 
        'Makefile': 1134
    }
    
    # We need a source of new timestamps for the 'touch' commands.
    time_counter = 1135
    def get_new_time():
        nonlocal time_counter
        t = time_counter
        time_counter += 1
        return t

    # 2. Define a function for the 'touch' command's action.
    def touch(filename):
        files[filename] = get_new_time()

    # 3. Represent the Makefile rules as a data structure.
    # The action is a function that will be called if the rule is triggered.
    rules = collections.OrderedDict([
        ('all',  {'deps': ['T', 'Z', 'X', 'Opps'], 'action': lambda: None}),
        ('T',    {'deps': ['Opps', 'X'], 'action': lambda: touch('A')}),
        ('X',    {'deps': ['Y'], 'action': lambda: touch('B')}),
        ('Z',    {'deps': ['Y'], 'action': lambda: touch('C')}),
        ('Opps', {'deps': ['T', 'Z'], 'action': lambda: touch('T')})
    ])

    # 4. Simulation engine
    already_built = set()
    build_stack = set()

    def build_target(target):
        # If we have already evaluated this target in this run, do not repeat.
        if target in already_built:
            return

        # Check for circular dependencies.
        if target in build_stack:
            # Real 'make' would print a warning here and drop the dependency. We just return.
            return
            
        build_stack.add(target)
        
        rule = rules.get(target)
        # If there is no rule, it is a source file (like 'Y') or an existing file.
        if not rule:
            build_stack.remove(target)
            already_built.add(target)
            return

        # Recursively build all dependencies first.
        for dep in rule['deps']:
            build_target(dep)

        # Determine if the current target needs to be rebuilt.
        needs_rebuild = False
        target_timestamp = files.get(target, -1)
        
        if target_timestamp == -1:
            # Rule runs if the target file does not exist.
            needs_rebuild = True
        else:
            # Rule runs if any dependency is newer than the target.
            for dep in rule['deps']:
                dep_timestamp = files.get(dep)
                if dep_timestamp is None:
                    # Also rebuild if a dependency file doesn't exist.
                    needs_rebuild = True
                    break
                if dep_timestamp > target_timestamp:
                    needs_rebuild = True
                    break

        # Execute the rule's action if rebuilding is necessary.
        if needs_rebuild:
            rule['action']()

        build_stack.remove(target)
        already_built.add(target)

    # Start the simulation from the 'all' target's dependencies.
    for dependency in rules['all']['deps']:
        build_target(dependency)
        
    # The 'ls' command in the 'all' rule shows the final state of the directory.
    final_files = sorted(files.keys())
    print("Final files list:")
    print(final_files)

# Run the simulation
simulate_make_all()