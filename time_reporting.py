import os


class TimeReportsKTG:

    def __init__(self, workspace, options):

        self.times_dictionary = None
        self.start_time = None

        self.workspace = workspace
        self.options = options

    def get_start_time(self):

        start_time = 0
        for heading, time in self.times_dictionary.items():
            if 'start' in heading.lower():
                start_time = time
                break
        self.start_time = start_time

    def append_dictionary(self, times):

        if not self.times_dictionary:
            self.times_dictionary = times
        else:
            self.times_dictionary = self.times_dictionary.update(times)

    def time_reporter(self, times, new_iter):

        # Check for existing values
        self.append_dictionary(times)
        if not self.start_time:
            self.get_start_time()

        time_results = os.path.join(self.workspace, "time_results.txt")
        time_printout = open(time_results, "a")
        if new_iter:
            time_printout.writelines(f"\n\n    ---- THIS IS A NEW ITERATION ----\n\n")
            for option_def, option in self.options.items():
                time_printout.writelines(f'\n{option_def}: {option}')
        time_printout.close()

        timenames = list(times.keys())
        elapsed_times = {}
        self.times_dictionary = self.times_dictionary.update(times)

        for timename in timenames:

            total_seconds = times[timename]
            elapsed = total_seconds - self.start_time
            hours = elapsed // 3600
            elapsed = elapsed - 3600 * hours
            minutes = elapsed // 60
            elapsed = elapsed - 60 * minutes
            seconds = round(elapsed, 2)
            elapsed_times[timename] = (hours, minutes, seconds)
            print(f"{timename} processing finished after: {hours} hours, {minutes} minutes, {seconds} seconds")
            time_printout = open(time_results, "a")
            time_printout.writelines(f"\n{timename}: \n{hours} hours, {minutes} minutes, {seconds} seconds\n")
            time_printout.close()
