/* [[file:~/Documents/work/qmckl_myFork/qmckl/org/qmckl_distance.org::*Headers][Headers:2]] */
#include "qmckl.h"
#include "assert.h"
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
int main() {
  qmckl_context context;
  context = qmckl_context_create();
/* Headers:2 ends here */

/* [[file:~/Documents/work/qmckl_myFork/qmckl/org/qmckl_distance.org::*Test][Test:2]] */
qmckl_exit_code test_qmckl_distance_sq(qmckl_context context);
assert(test_qmckl_distance_sq(context) == QMCKL_SUCCESS);
/* Test:2 ends here */

/* [[file:~/Documents/work/qmckl_myFork/qmckl/org/qmckl_distance.org::*Test][Test:2]] */
qmckl_exit_code test_qmckl_dist(qmckl_context context);
assert(test_qmckl_dist(context) == QMCKL_SUCCESS);
/* Test:2 ends here */

/* [[file:~/Documents/work/qmckl_myFork/qmckl/org/qmckl_distance.org::*End of files][End of files:1]] */
assert (qmckl_context_destroy(context) == QMCKL_SUCCESS);

  return 0;
}
/* End of files:1 ends here */
