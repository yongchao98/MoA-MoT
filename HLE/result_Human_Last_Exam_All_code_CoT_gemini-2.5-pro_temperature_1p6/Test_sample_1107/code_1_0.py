def simulate_make_all():
    """
    This script simulates the execution of 'make all' to determine the final set of files.
    The logic follows the behavior of GNU Make, including handling dependencies,
    timestamps, non-existent targets, and circular dependencies.
    """
    # Initial files in the directory. We use a set to track file existence.
    files = {'X', 'Y', 'Z', 'OPPS', 'Makefile'}

    # 1. Processing 'T' as the first dependency of 'all'.
    # 1a. For 'T', process dependency 'Opps'.
    #    - The rule 'Opps: T Z' has a circular dependency, so 'T' is ignored as a prereq.
    #    - Dependency 'Z' is up-to-date relative to 'Y'.
    #    - Target file 'Opps' doesn't exist, so its recipe 'touch T' runs.
    files.add('T')

    # 1b. For 'T', process dependency 'X'.
    #    - 'X' is out-of-date relative to 'Y', so its recipe 'touch B' runs.
    files.add('B')

    # 1c. Target 'T' itself is rebuilt because its prerequisite 'X' was remade.
    #    - The recipe 'touch A' runs.
    files.add('A')

    # 2. Processing 'Z', the second dependency of 'all'.
    #    - 'Z' is up-to-date, so nothing happens.

    # 3. Processing 'X', the third dependency of 'all'.
    #    - The recipe for 'X' has already been run in this invocation, so it's skipped.

    # 4. Processing 'Opps', the fourth dependency of 'all'.
    #    - The target file 'Opps' still does not exist, so its recipe 'touch T' must run again.
    files.add('T') # Updates timestamp, file is already in the set.

    # 5. The recipe for 'all' is 'ls', which shows the final files.
    final_files = sorted(list(files))

    print("The final files in the directory are:")
    for file_name in final_files:
        print(file_name)

simulate_make_all()