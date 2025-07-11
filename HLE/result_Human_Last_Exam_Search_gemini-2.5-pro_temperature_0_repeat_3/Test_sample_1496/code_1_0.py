def beaufort_decrypt(ciphertext, key):
    alphabet = "abcdefghijklmnopqrstuvwxyz "
    char_to_int = {char: i for i, char in enumerate(alphabet)}
    int_to_char = {i: char for i, char in enumerate(alphabet)}
    n = len(alphabet)
    
    plaintext = ""
    key_index = 0
    
    for char in ciphertext:
        if char in char_to_int:
            c = char_to_int[char]
            k_char = key[key_index % len(key)]
            k = char_to_int[k_char]
            
            p = (k - c + n) % n
            plaintext += int_to_char[p]
            
            key_index += 1
        else:
            plaintext += char
            
    return plaintext

ciphertext = "fetssonayhmcippuqadksfd dbhmquecmomiaqqadmzk lfvjqmydbumamsaaodqdjxozpr fexhiiiqfgyxjoapdadaygq idmtapibftejjnvlanmksrbzfijteknbpslxehmlkuqddyeixsdfbinqxlf wgxuwaejssyzuzmlumitoumfwhpmzbokgfn wsvllbmbfoyushhglfly"
key = "williamtellwoollennightcapsimeoni"

decrypted_text = beaufort_decrypt(ciphertext, key)
print(decrypted_text)