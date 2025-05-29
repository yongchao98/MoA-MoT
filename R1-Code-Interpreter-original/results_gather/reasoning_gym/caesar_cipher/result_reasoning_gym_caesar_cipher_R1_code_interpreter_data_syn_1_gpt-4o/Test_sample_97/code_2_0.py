# The encrypted message with shift 4
shift_4_text = "GUR VAINYVQVGL BE HARASBEPRNOVYVGL BS NAL CEBIVFVBA BS GUVF NTERRZRAG FUNYY ABG IBVQ GUR ERZNVAVAT CEBIVFVBAF"

# Decrypt using ROT13 (shift of 13)
decrypted_message = decrypt_caesar(shift_4_text, 13)
print(decrypted_message)