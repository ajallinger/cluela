/*######################################################################
32 bit Fowler/Noll/Vo FNV-1a hash code
Public Domain
########################################################################
changed return type to INTEGER*4            2017/04/04
translated to FORTRAN                       2013, Andy Allinger
                                            andy_a@users.sourceforge.net

Revision: 5.1                               2009/06/30 09:13:32
########################################################################
Fowler/Noll/Vo hash

The basis of this hash algorithm was taken from an idea sent
as reviewer comments to the IEEE POSIX P1003.2 committee by:

     Phong Vo (http://www.research.att.com/info/kpv/)
     Glenn Fowler (http://www.research.att.com/~gsf/)

 In a subsequent ballot round:

     Landon Curt Noll (http://www.isthe.com/chongo/)

improved on their algorithm.  Some people tried this hash
and found that it worked rather well.  In an EMail message
to Landon, they named it the ``Fowler/Noll/Vo'' or FNV hash.

FNV hashes are designed to be fast while maintaining a low
collision rate. The FNV speed allows one to quickly hash lots
of data while maintaining a reasonable collision rate.  See:

     http://www.isthe.com/chongo/tech/comp/fnv/index.html

for more details as well as other forms of the FNV hash.

#######################################################################
Please do not copyright this code.  This code is in the public domain.

LANDON CURT NOLL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
EVENT SHALL LANDON CURT NOLL BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

By:
	chongo <Landon Curt Noll> /\oo/\
     http://www.isthe.com/chongo/

Share and Enjoy!	:-)

#######################################################################
32 bit FNV-1 and FNV-1a non-zero initial basis

The FNV-1 initial basis is the FNV-0 hash of the following 32 octets:

             chongo <Landon Curt Noll> /\../\

NOTE: The \'s above are not back-slashing escape characters.
They are literal ASCII  backslash 0x5c characters.

NOTE: The FNV-1a initial basis is the same value as FNV-1 by definition.

#######################################################################
perform a 32 bit Fowler/Noll/Vo FNV-1a hash on a buffer

input:
	BUFFER	- start of buffer to hash
	LENGTH	- length of buffer in octets [bytes]
	HASH  	- previous hash value or [standard initial value]
 assumption:  INTEGER's are at least 32 bits

NOTE:  To use the recommended 32 bit FNV-1a hash, use the initial value:
              HASH = -2128831035
       as the argument on the first call to FNV32

!NOTE:  Only pass the previous hash value if you always hash the same
        data in the same order (or else you will get different answers!)

output:
	32 bit HASH

######################################################################*/
#include <stdint.h>
void fnv32 (void *buffer, int length, int32_t *hash)
{
	unsigned char *bp = (unsigned char *) buffer;
	int j;

/* FNV-1a hash each octet in the buffer */
	for (j = 0; j < length; ++j) {
		*hash ^= bp[j];
		*hash *= 16777619;  /* multiply by the 32 bit FNV magic prime */
	}
	return;
} /* End of fnv32 */
