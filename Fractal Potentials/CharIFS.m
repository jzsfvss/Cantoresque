% Character of an IFS.

function chr = CharIFS(phi, prb)

chr = sum(prb.*log(phi));
