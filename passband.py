# %%
import phoebe
from phoebe import u # units

# %%
ptf = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/SDSS.g.ptf.txt'
pb = phoebe.atmospheres.passbands.Passband(ptf=ptf, pbset='SDSS', pbname='g', effwl=4671.78, wlunits=u.AA, calibrated=True, reference='SVO Filter Profile Service', version=1.0, comments='')

# %%
pb.compute_blackbody_response()
pb.compute_bb_reddening(verbose=True)
print(pb.content)

# %%
# ck2004 = 

pb.compute_ck2004_response(path=ck2004, verbose=True)
pb.compute_ck2004_intensities(path=ck2004, verbose=True)
pb.compute_ck2004_ldcoeffs()
pb.compute_ck2004_ldints()
pb.compute_ck2004_reddening(path=ck2004, verbose=True)

# pb.compute_phoenix_response(path='tables/phoenix', verbose=True)
# pb.compute_phoenix_intensities(path='tables/phoenix', verbose=True)
# pb.compute_phoenix_ldcoeffs()
# pb.compute_phoenix_ldints()
# pb.compute_phoenix_reddening(path='tables/phoenix', verbose=True)

pb.save('SDSS.g.pb')
# %%
