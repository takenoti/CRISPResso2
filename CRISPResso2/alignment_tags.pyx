from cpython.unicode cimport PyUnicode_AsUTF8, PyUnicode_FromStringAndSize

def aligned_crispresso_to_cs_tag(str reference, str read, bint short_cs_tag=False):
    cdef int i = 0
    cdef int sequence_len = len(reference)
    cdef int match_start, insertion_start, deletion_start
    cdef char dash = b'-'  # dash character in C
    cdef list operations = []

    # Obtain C pointers to the UTF-8 encoded data of the Unicode strings.
    cdef const char* reference_ptr = PyUnicode_AsUTF8(reference)
    cdef const char* read_ptr = PyUnicode_AsUTF8(read)

    cdef char substitution[3]
    substitution[0] = b'*'  # first character for a substitution operation

    cdef char reference_char, read_char  # local variables for current characters

    while i < sequence_len:
        reference_char = reference_ptr[i]
        read_char = read_ptr[i]

        # If characters match, fast-forward over all matching characters.
        if reference_char == read_char:
            match_start = i
            i += 1
            while i < sequence_len and reference_ptr[i] == read_ptr[i]:
                i += 1

            if short_cs_tag:
                # Append a match operation: ":" followed by the matching len.
                operations.append(":" + str(i - match_start))
            else:
                # Append a match operation: "=" followed by the matching substring.
                operations.append("=" + reference[match_start:i])
            continue

        # Insertion: the reference has a dash.
        if reference_char == dash:
            insertion_start = i
            i += 1
            while i < sequence_len and reference_ptr[i] == dash:
                i += 1
            operations.append("+" + read[insertion_start:i])
            continue

        # Deletion: the read has a dash.
        if read_char == dash:
            deletion_start = i
            i += 1
            while i < sequence_len and read_ptr[i] == dash:
                i += 1
            operations.append("-" + reference[deletion_start:i])
            continue

        # Substitution: use the '*' operation followed by the two differing characters.
        substitution[1] = reference_char
        substitution[2] = read_char
        operations.append(PyUnicode_FromStringAndSize(substitution, 3))
        i += 1

    return "".join(operations)
