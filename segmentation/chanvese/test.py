from chanvese import chanvese, mask2phi
import json
import numpy as np
from PIL import Image
from skimage.segmentation import chan_vese
import scipy.ndimage as nd

if __name__ == "__main__":
    with open('/Users/tinnguyen/Work/pattern_selection/speed.json', 'r') as f:
        speed = np.array(json.load(f))

    un = speed / np.max(speed)
    im = Image.fromarray(np.uint8(un*255), 'L')
    im.save('speed.jpg')

    # init_mask = speed < 50
    # # init_phi = mask2phi(init_mask)
    # init_phi = nd.distance_transform_edt(speed > 50)
    # un = speed / np.max(speed)
    # congestion_mask = chan_vese(un, mu=0, lambda1=1, lambda2=1, tol=1e-3,
    #                         max_iter=500,dt=0.5, init_level_set=init_phi, extended_output=True)
    # im = Image.fromarray(np.uint8(congestion_mask[0])*255, 'L')
    # im.save('congestion_mask_skimage.jpg')



    # seg, phi, its = chanvese(speed, init_mask=speed<50)

    # im = Image.fromarray(np.uint8(seg)*255, 'L')
    # im.save('congestion_mask.jpg')