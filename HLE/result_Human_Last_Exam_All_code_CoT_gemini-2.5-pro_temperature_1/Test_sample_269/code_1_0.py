def cosmo_loco_output_stream():
    """
    This generator models the output of the Scheme program.
    The program first prints 1, then 1, then 2, and then
    enters an infinite loop printing 2.
    """
    yield 1
    yield 1
    yield 2
    while True:
        yield 2

# Let's print the first 15 numbers from the stream to see the pattern.
output_generator = cosmo_loco_output_stream()
output_list = [str(next(output_generator)) for _ in range(15)]
print(" ".join(output_list) + " ...")
