def optimized_dismantle_and_combine():
    # Initial counts
    machines_A = 3
    machines_B = 2
    machines_C = 2
    parts_X = 0
    parts_Y = 0
    parts_Z = 0

    # Convert Machines B and C into Machine A
    convert_to_A = min(machines_B, machines_C)
    machines_A += convert_to_A
    machines_B -= convert_to_A
    machines_C -= convert_to_A

    # Use Machines A to produce Parts X and Y
    parts_X += 2 * machines_A
    parts_Y += machines_A
    machines_A = 0

    # Use remaining Machines B to produce Part X
    parts_X += machines_B // 2
    machines_B %= 2

    # Use remaining Machines C to produce Part Y
    parts_Y += machines_C // 2
    machines_C %= 2

    # Convert Parts X and Y into Part Z
    convert_to_Z = min(parts_X, parts_Y)
    parts_Z += convert_to_Z
    parts_X -= convert_to_Z
    parts_Y -= convert_to_Z

    # Output the final counts
    result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
    print(result)

optimized_dismantle_and_combine()