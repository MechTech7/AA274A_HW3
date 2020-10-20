#!/usr/bin/env python

import numpy as np
import time
import cv2
import matplotlib.pyplot as plt

'''Error Output
('Iteration: ', 1)
('Vert padding: ', 1)
('Horiz padding: ', 100)
(579, 968, 3)
('F_shape: ', (3, 200, 3))
('G_shape: ', (576, 768))
(3, 201, 3)
Traceback (most recent call last):
  File "linear_filter.py", line 155, in <module>
    main()
  File "linear_filter.py", line 149, in main
    corr_img = corr(filt, test_card)
  File "linear_filter.py", line 49, in corr
    pix_val = np.dot (img_val, flat_filter)
ValueError: shapes (1809,) and (1800,) not aligned: 1809 (dim 0) != 1800 (dim 0)
'''
def corr(F, I):
    """
    Input
        F: A (k, ell, c)-shaped ndarray containing the k x ell filter (with c channels).
        I: An (m, n, c)-shaped ndarray containing the m x n image (with c channels).

    Returns
        G: An (m, n)-shaped ndarray containing the correlation of the filter with the image.
    """
    ########## Code starts here ##########
    
    k = F.shape[0]
    l = F.shape[1]

    k_f = k // 2
    l_f = l // 2

    k_r = k - k_f
    l_r = l - l_f

    padded_input = np.pad(I, pad_width=((k_f, k_r), (l_f, l_r), (0, 0)), mode='constant')
    flat_filter = np.ravel(F)

    print ("Vert padding: ", k_f)
    print ("Horiz padding: ", l_f)

    #print (padded_input)
    print(padded_input.shape)
    
    G = np.zeros((I.shape[0], I.shape[1]))
    print ("F_shape: ", F.shape)
    print ("G_shape: ", G.shape)
    for i in range(I.shape[0]):
        for j in range(I.shape[1]):
            p_i = i + k_f
            p_j = j + l_f

            image_region = padded_input[p_i - k_f:p_i + k_f+1, p_j - l_f:p_j + l_f+1]
            print (image_region.shape)
            
            print (image_region[0, 0])
            print (image_region[0, -1])
            img_val = np.ravel(image_region)
            pix_val = np.dot (img_val, flat_filter)
            G[i, j] = pix_val
    
    return G
    #raise NotImplementedError("Implement me!")
    ########## Code ends here ##########


def norm_cross_corr(F, I):
    """
    Input
        F: A (k, ell, c)-shaped ndarray containing the k x ell filter (with c channels).
        I: An (m, n, c)-shaped ndarray containing the m x n image (with c channels).

    Returns
        G: An (m, n)-shaped ndarray containing the normalized cross-correlation of the filter with the image.
    """
    ########## Code starts here ##########
    k = F.shape[0]
    l = F.shape[1]

    k_f = k // 2
    l_f = l // 2

    k_r = k - k_f
    l_r = l - l_f

    padded_input = np.pad(I, pad_width=((k_f, k_r), (l_f, l_r), (0, 0)), mode='constant')
    #padded_input = np.pad(I, pad_width=((k_f, k_r), (l_f, l_r)), mode='constant')
    flat_filter = np.ravel(F)

    filt_norm = np.linalg.norm(flat_filter)
    #print (padded_input)
    print(padded_input.shape)
    
    G = np.zeros_like((I.shape[0], I.shape[1]))
    for i in range(I.shape[0]):
        for j in range(I.shape[1]):
            p_i = i + k_f
            p_j = j + l_f

            img_val = np.ravel(padded_input[p_i - k_f:p_i + k_f+1, p_j - l_f:p_j + l_f+1])
            pix_val = np.dot (img_val, flat_filter)
            
            im_norm = np.linalg.norm(img_val)

            pix_val /= (im_norm * filt_norm)
            G[i, j] = pix_val


    return G

    #raise NotImplementedError("Implement me!")
    ########## Code ends here ##########


def show_save_corr_img(filename, image, template):
    # Not super simple, because need to normalize image scale properly.
    fig, ax = plt.subplots()
    cropped_img = image[:-template.shape[0], :-template.shape[1]]
    im = ax.imshow(image, interpolation='none', vmin=cropped_img.min())
    fig.colorbar(im)
    fig.savefig(filename, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def main():
    test_card = cv2.imread('test_card.png').astype(np.float32)

    image = np.array([[[1], [2], [3]], [[4], [5], [6]], [[7], [8], [9]]], dtype=np.float32)
    
    filt1 = np.zeros((3, 3, 1))
    filt1[1, 1] = 1

    filt2 = np.zeros((3, 200, 1))
    filt2[1, -1] = 1

    filt3 = np.zeros((3, 3, 1))
    filt3[:, 0] = -1
    filt3[:, 2] = 1

    filt4 = (1./273.)*np.array([[1, 4, 7, 4, 1],
                              [4, 16, 26, 16, 4],
                              [7, 26, 41, 26, 7],
                              [4, 16, 26, 16, 4],
                              [1, 4, 7, 4, 1]])
    filt4 = np.expand_dims(filt4, -1)

    grayscale_filters = [filt1, filt2, filt3, filt4]

    color_filters = list()
    for filt in grayscale_filters:
        # Making color filters by replicating the existing
        # filter per color channel.
        color_filters.append(np.concatenate([filt, filt, filt], axis=-1))

    for idx, filt in enumerate(color_filters):
        print ("Iteration: ", idx)
        start = time.time()
        corr_img = corr(filt, test_card)
        stop = time.time()
        print 'Correlation function runtime:', stop - start, 's'
        show_save_corr_img("corr_img_filt%d.png" % idx, corr_img, filt)

if __name__ == "__main__":
    main()
