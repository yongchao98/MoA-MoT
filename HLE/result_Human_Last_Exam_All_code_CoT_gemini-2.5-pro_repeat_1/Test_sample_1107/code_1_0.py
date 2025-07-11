import sys

def solve():
    """
    This script simulates the behavior of the 'make' command on the given Makefile
    to determine the final set of files in the directory.
    """
    # Timestamps are represented as integers. Higher is newer.
    files = {
        'X': 1051,
        'Y': 1052,
        'Z': 1054,
        'OPPS': 1132,
        'Makefile': 1134,
    }

    # Makefile rules
    rules = {
        'all': {'deps': ['T', 'Z', 'X', 'Opps'], 'command': 'ls'},
        'T': {'deps': ['Opps', 'X'], 'command': 'touch A'},
        'X': {'deps': ['Y'], 'command': 'touch B'},
        'Z': {'deps': ['Y'], 'command': 'touch C'},
        'Opps': {'deps': ['T', 'Z'], 'command': 'touch T'},
        'Q': {'deps': ['T'], 'command': 'touch H'},
    }

    time_counter = 1135
    remade_this_run = set()

    def make(target, call_stack=None):
        nonlocal time_counter
        if call_stack is None:
            call_stack = set()

        if target in remade_this_run:
            return

        if target in call_stack:
            return  # Circular dependency detected

        call_stack.add(target)

        rule = rules.get(target)
        if not rule:
            if target not in files:
                # This would be an error in a real make execution
                print(f"Error: No rule to make target '{target}'", file=sys.stderr)
            call_stack.remove(target)
            return

        for dep in rule['deps']:
            make(dep, call_stack)

        must_rebuild = False
        if target not in files:
            must_rebuild = True
        else:
            target_time = files.get(target, 0)
            for dep in rule['deps']:
                if dep in remade_this_run:
                    must_rebuild = True
                    break
                dep_time = files.get(dep, 0)
                if dep_time > target_time:
                    must_rebuild = True
                    break
        
        if must_rebuild:
            command = rule['command']
            if command.startswith('touch'):
                file_to_touch = command.split(' ')[1]
                files[file_to_touch] = time_counter
                time_counter += 1
            remade_this_run.add(target)

        call_stack.remove(target)

    # Start the simulation from the 'all' target
    make('all')

    # Print the sorted list of final files
    final_files = sorted(list(files.keys()))
    for f in final_files:
        print(f)

solve()