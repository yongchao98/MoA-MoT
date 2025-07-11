# First, you may need to install the Biopython library. You can do this by running:
# pip install biopython

import sys
import io
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: Biopython library not found.")
    print("Please install it using: pip install biopython")
    sys.exit(1)

# The amino acid sequence to identify
amino_acid_sequence = (
    "MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP "
    "TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED "
    "FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV "
    "TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH "
    "YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD "
    "EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN "
    "NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG "
    "VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY "
    "TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH "
    "FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF "
    "QYEGTYKWVN PHKL"
)

# Clean the sequence by removing whitespace
cleaned_sequence = "".join(amino_acid_sequence.split())

print("Performing BLAST search against NCBI's non-redundant protein database.")
print("This can take a minute, please be patient...")

final_answer = ""
try:
    # Use qblast to send the sequence to the NCBI BLAST server
    # 'blastp' is for protein-protein comparison
    # 'nr' is the non-redundant protein database
    result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)

    # Read the result from the handle and treat it as a string
    blast_xml_string = result_handle.read()
    result_handle.close()
    
    # To parse the XML, we can treat the string as a file using io.StringIO
    blast_in = io.StringIO(blast_xml_string)

    # Check if any result was returned
    if not blast_xml_string.strip():
         final_answer = "Search returned no result. The sequence might be invalid or too short."
    else:
        # Parse the BLAST result using NCBIXML
        blast_record = NCBIXML.read(blast_in)

        # Check if the search yielded any hits (alignments)
        if len(blast_record.alignments) > 0:
            # The first alignment is the best match
            best_alignment = blast_record.alignments[0]
            protein_title = best_alignment.title
            
            print(f"\n--- Search Complete ---\n")
            print(f"Top BLAST Hit Full Title: {protein_title}\n")

            # Try to extract the most descriptive name, often found after "RecName: Full="
            if "RecName: Full=" in protein_title:
                # e.g., RecName: Full=Dihydrolipoyl dehydrogenase, mitochondrial; Short=DLD
                final_answer = protein_title.split("RecName: Full=")[1].split(';')[0]
            elif '[' in protein_title:
                # Fallback to taking everything before the organism name
                # e.g., dihydrolipoamide dehydrogenase [Homo sapiens]
                final_answer = protein_title.split('[')[0].strip().split(' ', 1)[1]
            else:
                # Final fallback to just the title (without the accession code)
                final_answer = protein_title.split(' ', 1)[1] if ' ' in protein_title else protein_title
        else:
            final_answer = "No matching protein found."

except Exception as e:
    final_answer = f"An error occurred: {e}\nPlease check your internet connection."

print(f"Identified Protein: {final_answer.strip()}")