import matplotlib.pyplot as plt
import cv2


def main():
    image = cv2.imread('helix.jpg')
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    blue_hue = 107
    blue_mask = cv2.inRange(hsv, (98, 70, 0), (115, 255, 255))
    red_hue = 0
    red1_mask = cv2.inRange(hsv, (0, 120, 0), (9, 255, 255))
    red2_mask = cv2.inRange(hsv, (150, 30, 0), (180, 255, 255))
    red_mask = red1_mask | red2_mask
    # red_mask = cv2.inRange(hsv, (1, 0, 0), (0, 0, 0))
    orange_mask = cv2.inRange(hsv, (110, 60, 20), (16, 255, 220))
    blue_masked = cv2.bitwise_and(hsv, hsv, mask=blue_mask)
    red_masked = cv2.bitwise_and(hsv, hsv, mask=red_mask)
    orange_masked = cv2.bitwise_and(hsv, hsv, mask=orange_mask)
    bg_masked = cv2.bitwise_and(hsv, hsv, mask=~(blue_mask |
                                                 red_mask |
                                                 orange_mask))

    h, s, v = cv2.split(blue_masked)
    hue_shift = blue_hue - red_hue
    # hue_shift *= 0
    shifting_down = h > hue_shift
    h[shifting_down] -= hue_shift
    h[~shifting_down] += 180 - hue_shift
    blue_masked = cv2.merge((h, s, v))
    blue_masked = cv2.bitwise_and(blue_masked, blue_masked, mask=blue_mask)

    h, s, v = cv2.split(red_masked)
    shifting_up = h < 180 - hue_shift
    h[shifting_up] += hue_shift
    h[~shifting_up] -= 180 - hue_shift
    red_masked = cv2.merge((h, s, v))
    red_masked = cv2.bitwise_and(red_masked, red_masked, mask=red_mask)
    # red_masked *= 0
    # blue_masked *= 0
    # orange_masked *= 0

    masked = red_masked | blue_masked | orange_masked | bg_masked
    masked = cv2.flip(masked, 1)
    masked = cv2.cvtColor(masked, cv2.COLOR_HSV2BGR)
    cv2.imwrite('swapped.jpg', masked)
    masked = cv2.cvtColor(masked, cv2.COLOR_BGR2RGB)

    plt.imshow(masked)
    plt.gca().set_axis_off()
    plt.tight_layout(0)
    plt.show()


main()
