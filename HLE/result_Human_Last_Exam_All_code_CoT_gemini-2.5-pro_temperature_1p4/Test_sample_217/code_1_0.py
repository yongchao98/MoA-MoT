def find_correct_triplets():
    """
    This function identifies and prints the three correct triplets from the user's list.
    A triplet is considered correct if the two species listed are both clear examples
    of the associated mode of mimicry, crypsis, or warning signaling.
    """
    
    # The list of all triplets provided by the user.
    triplets = {
        1: "1) Aristotelian, Charadrius vociferus - wing, Recurvirostra avosetta - wing",
        2: "2) Automimicry, Arcas cypria - tail, Apis melifera - abdomen",
        3: "3) Batesian, Eristalini - color, Melipotini - color",
        4: "4) Gilbertian, Passiflora - leaf, Myrmecia chrysogaster - venom",
        5: "5) Müllerian, Heliconiini - color, Melipotini - color",
        6: "6) Vavilovian, Secale cereale - seeds, Asclepias speciosa - seeds",
        7: "7) Camouflage, Liturgusa maya - abdomen, Limenitis archippus - wing",
        8: "8) Aposematism, Danaus plexipus - wing, Cycnia tenera - tymbal"
    }
    
    # The numbers of the identified correct triplets, in ascending order.
    correct_indices = [1, 5, 8]
    
    print("The three correct triplets, in ascending order, are:")
    for index in correct_indices:
        print(f"\n{triplets[index]}")

if __name__ == '__main__':
    find_correct_triplets()