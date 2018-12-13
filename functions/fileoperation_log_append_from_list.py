###This is a very simple program that writes out a generic log, with multiple
###levels indicated by the number of #s at the start of a line. This may be
###useful for troubleshooting or for saving simple results.

###GSJ 27/09/2017

def fileoperation_log_append_from_list(logfile, line_list = [[]]):
    f = open(logfile, 'ab')
    for indent in range(len(line_list)):
        for line in line_list[indent]:
            line_to_write = ['#' for i in range(indent)] + [line] + ['\n']
            line_to_write = ''.join(line_to_write)
            f.write(line_to_write)
    f.close()
    return None
