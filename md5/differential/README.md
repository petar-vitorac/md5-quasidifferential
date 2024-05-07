# MD5 - Differential

- `md5coll.c` finds collisions in MD5 that follow Wang's characteristic. It has to be compiled for 32-bit (`-m32` flag).
- `md5.sage` simulates the execution of the MD5 function. During the execution, it can also extract intermediary differences, or (depending on the mode) check whether the message follows Wang's characteristic. The first output of `md5coll.c` can be inserted in `md5.sage` as `M0` to perform the calculations on that message.
