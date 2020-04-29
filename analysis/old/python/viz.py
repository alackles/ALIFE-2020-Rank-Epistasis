import pygame
import time

def visualize(sorted_ids, mut_sorted_ids, delay = 1, sw = 800, sh = 800, tw = 7, th = 7):
    screen_width = sw 
    screen_height = sh 
    tile_width = tw 
    tile_height = th    
    
    pygame.init()
    screen = pygame.display.set_mode((screen_width, screen_height))

    screen.fill((0,0,0))
    for i in range(len(sorted_ids)):
        old_rank = i
        org_id = sorted_ids[i]
        new_rank = 0
        for j in range(len(sorted_ids)):
            if mut_sorted_ids[j] == org_id:
                new_rank = j
                break
        color = (150,150,150)
        if new_rank < old_rank:
            color = (0,150,0)
        elif new_rank > old_rank:
            color = (150,0,0)
        pygame.draw.line(\
            screen, \
            color, \
            (100,                old_rank * (tile_height + 1) + (tile_height / 2)), \
            (screen_width - 100, new_rank * (tile_height + 1) + (tile_height / 2)),
            1)   
    for i in range(len(sorted_ids)):
        pygame.draw.rect(screen, \
            (150,150,150), \
            (100, i * (tile_height + 1), tile_width, tile_height))
        pygame.draw.rect(screen, \
            (150,150,150), \
            (screen_width - 100, i * (tile_height + 1), tile_width, tile_height))

    pygame.display.flip()
    time.sleep(delay)
    pygame.quit()
