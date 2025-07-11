import sys

def solve():
    """
    This function analyzes the two questions about the toy instruction set and prints the answers.
    """

    # --- Analysis for Question A ---
    # Question A asks if every boopable function is Turing-computable.
    # A function is "boopable" if there exists a program L that can compute it
    # for any input x, given that the machine's resources (R, A, W) can be
    # made arbitrarily large.
    #
    # A machine model with an instruction set featuring arithmetic (ADD),
    # bitwise logic (BITWISENAND), memory access (LOAD/STORE), and conditional
    # branching (BRANCHIFZERO), coupled with unbounded memory and registers,
    # is known to be Turing-complete. It is equivalent to a Random Access Machine (RAM).
    #
    # Since the machine is Turing-complete, any function it can compute
    # (i.e., any "boopable" function) must be, by definition, a Turing-computable function.
    answer_A = "Yes"

    # --- Analysis for Question B ---
    # Question B asks if a specific program L can exist under fixed constraints:
    # (R,A,W)=(64,64,512), len(L) <= 2^32+50, input x in [0, 2^32),
    # output x boops, and total steps <= x + 50.
    #
    # The time constraint `steps <= x + 50` is the key. For large x, this
    # implies that each boop must be generated in approximately one machine step.
    # This rules out any simple iterative loop, as a loop requires instructions
    # for the operation (BOOP), counter update (ADD), and condition checking/jumping
    # (BRANCHIFZERO), making each iteration take more than one step.
    #
    # The only way to achieve one step per boop is to have a long, contiguous
    # block of `BOOP` instructions in the program (an "unrolled loop"). To get
    # exactly x boops, the program must jump to a specific starting instruction
    # within this block, based on a calculation involving x. This is a "computed jump".
    #
    # The instruction set's branch instruction is `BRANCHIFZERO adr, reg`, where
    # the destination `adr` is a fixed immediate value, not a value from a register.
    # The architecture lacks a "jump to register" instruction.
    #
    # Furthermore, it's not possible to efficiently simulate a computed jump.
    # Self-modification is infeasible because the `STORE adr <- reg` instruction also
    # uses a fixed immediate address. A dispatch tree of branches would take O(x)
    # steps, violating the time constraint.
    #
    # Since there is no way to enter an unrolled loop at a calculated offset
    # within the given time budget, no such program L can exist.
    answer_B = "No"

    print(f"{answer_A} {answer_B}")

solve()