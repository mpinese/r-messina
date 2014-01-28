/* RInterrupt.h
 * A hack to support user interruption of C++ code from R.  Code is by Simon
 * Urbanek, from https://stat.ethz.ch/pipermail/r-devel/2011-April/060702.html.
 * 
 */

namespace RCheckInterrupt
{

static void chkIntFn(void *dummy)
{
	R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
bool checkInterrupt() 
{
	return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

}
