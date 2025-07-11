import sys

def solve():
  """
  This function encapsulates the reasoning for the two-part question.

  Part A: Is every boopable function Turing-computable?
  A boopable function is defined by a program `L` that works on a machine
  where the parameters R (registers), A (address space), and W (word size)
  can be made arbitrarily large. A machine model with unbounded memory,
  combined with arithmetic, universal logic (NAND), and conditional branching,
  is Turing-complete. Therefore, any function computable by this machine model
  (a boopable function) must be Turing-computable. The answer is Yes.

  Part B: Does a specific fast program exist?
  The constraints are: machine (64,64,512), input x in [0, 2^32),
  program length <= 2^32+50, and time <= x+50.

  The time constraint `x+50` strongly implies a non-looping execution for the
  boops, specifically a long sequence of `BOOP` instructions (a "BOOP sled").
  The program must dynamically jump to the correct starting point in this sled
  to execute exactly `x` boops.

  This can be achieved with self-modifying code:
  1. The program L is structured with a setup section followed by a BOOP sled of
     length 2^32. The HALT instruction is implicitly at the end (memory is zero-initialized).
  2. The setup code calculates the target address to jump to:
     `target_addr = (sled_end_address) - x`.
  3. It then constructs a `BRANCHIFZERO target_addr, r_zero` instruction.
  4. Crucially, we assume a "favorable" instruction encoding (not ruled out by the
     problem description) where the address part of the instruction can be set
     with a simple addition, avoiding a slow bit-shift operation.
  5. The setup code then `STORE`s this instruction in memory and executes it.

  The setup time is a small constant (e.g., ~10 steps). The total time is
  `T_setup + T_jump + T_boops + T_halt`, which is approximately `10 + 1 + x + 1 = x + 12`.
  This is less than `x + 50`. Thus, such a program can exist. The answer is Yes.
  """
  answer_A = "Yes"
  answer_B = "Yes"
  print(f"{answer_A} {answer_B}")

solve()
<<<Yes Yes>>>