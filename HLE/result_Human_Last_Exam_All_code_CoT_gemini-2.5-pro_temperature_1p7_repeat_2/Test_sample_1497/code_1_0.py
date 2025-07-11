# This script prints the solution based on the analysis of how different playing variations
# affect the physical parameters of a violin's sound production.

# (1) sul ponticello -> changes bowing, affecting string harmonic damping -> group ii
# (2) bridge mute -> changes body resonance, especially higher frequencies -> group iv
#     - a mute adds mass, lowering the resonance frequency f2 -> down
# (3) helium environment -> changes air resonance frequency -> group iii
# (4) changing from A string to E string -> changes fundamental frequency -> group i

answer = "ii,iv,iii,i,down"
print(answer)