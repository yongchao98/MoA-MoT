def beaufort_decrypt(ciphertext, key):
    alphabet = "abcdefghijklmnopqrstuvwxyz "
    decrypted_text = ""
    key_index = 0
    for char in ciphertext:
        if char in alphabet:
            c = alphabet.index(char)
            k = alphabet.index(key[key_index % len(key)])
            p = (k - c + len(alphabet)) % len(alphabet)
            decrypted_text += alphabet[p]
            key_index += 1
        else:
            decrypted_text += char
    return decrypted_text

key = "williamtellwoollennightcapsimeoni"
ciphertext = "fetssonayhmcippuqadksfd dbhmquecmomiaqqadmzk lfvjqmydbumamsaaodqdjxozpr fexhiiiqfgyxjoapdadaygq idmtapibftejjnvlanmksrbzfijteknbpslxehmlkuqddyeixsdfbinqxlf wgxuwaejssyzuzmlumitoumfwhpmzbokgfn wsvllbmbfoyushhglfly"

decrypted_message = beaufort_decrypt(ciphertext, key)
print(decrypted_message)
