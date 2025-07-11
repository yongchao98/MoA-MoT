def beaufort_decipher(ciphertext, key):
    charset = "abcdefghijklmnopqrstuvwxyz "
    deciphered_text = ""
    key_index = 0
    for char_c in ciphertext:
        if char_c in charset:
            c = charset.find(char_c)
            k = charset.find(key[key_index % len(key)])
            p = (k - c + len(charset)) % len(charset)
            deciphered_text += charset[p]
            key_index += 1
        else:
            deciphered_text += char_c
    return deciphered_text

ciphertext = "fetssonayhmcippuqadksfd dbhmquecmomiaqqadmzk lfvjqmydbumamsaaodqdjxozpr fexhiiiqfgyxjoapdadaygq idmtapibftejjnvlanmksrbzfijteknbpslxehmlkuqddyeixsdfbinqxlf wgxuwaejssyzuzmlumitoumfwhpmzbokgfn wsvllbmbfoyushhglfly"
key = "williamtellnightcapsimeoni"

deciphered_text = beaufort_decipher(ciphertext, key)
print(deciphered_text)