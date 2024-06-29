package bdmmflow.utils;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * This class is a LRU cache that can be used to memoize values. It only stores
 * cacheSize number of items and discards the oldest values in the cache.
 */
public class LRUCache<K, V> extends LinkedHashMap<K, V> {
    private final int cacheSize;

    public LRUCache(int cacheSize) {
        super(20, 0.75F, false);
        this.cacheSize = cacheSize;
    }

    @Override
    protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
        return size() >= cacheSize;
    }
}
