package bdmmflow.utils;

import java.util.LinkedHashMap;
import java.util.Map;

public class LRUCache<K, V> extends LinkedHashMap<K, V> {
    private final int cacheSize;

    public LRUCache(int cacheSize) {
        super(20, 0.75F, false);
        this.cacheSize = cacheSize;
    }

//    protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
//        return size() >= cacheSize;
//    }
}
