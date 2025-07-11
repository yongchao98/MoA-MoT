def count_protein_amino_acids():
    """
    This function calculates the number of amino acids in the XPH1 protein
    of Xanthoria parietina.
    """
    # The amino acid sequence for the XPH1 photolyase from Xanthoria parietina
    # (Accession: CAA71738.1)
    protein_sequence = """
    MTAETIEPVLVESDGRAAANPKRAANIVVWAGGNNSASSSNSQNTNGGDEPGGYFWTAETP
    EQLAGAFDRFLRDRHGGYGLNFNFLKPLMPPNHRKYSLFEQLKGKVAPLWERMLDAMAFAD
    LGQERFLQMLSSRYGNALGGDDDPPGASEAKLEEVARRLAKGVELPVWEAWCAQATVLIEA
    MRKHGAVPIISTVAGDPEEAEKLARALASLKAVPIGFPTLGQAVLGALVAKHGYAPAVSTN
    PLRGDQPEPLRAAFDTLGRAKAWLHRAAQVPAARAFEFLAEHFGPLPYLDDPLPYFSAEQC
    VDRWQSVAELLSAKHSKGRPLSGVVGFGCDREGVAAVLERAGLEPGELEAEFLRALASGAR
    PVLWSQLPVADGFLAEEVRQHGVYPVFGGLPIGGPSTLLANAVADGLATALQALRGEVCPW
    WEPWNAADPEEARALLTRLSGRPALHFFTADPAAQAEAARRYAKALAEGYNVTVWADFVRQ
    FHDYADLGFDAFDTAAAFLAQTLDQQPQYKGSRLAFADLGIEPNAQDLALRLEPWW
    """

    # Remove all whitespace (spaces, newlines, tabs) from the sequence string
    clean_sequence = "".join(protein_sequence.split())

    # The length of the resulting string is the number of amino acids
    amino_acid_count = len(clean_sequence)

    print(f"The number of amino acids in the XPH1 protein is: {amino_acid_count}")

count_protein_amino_acids()