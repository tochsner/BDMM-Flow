package bdmmflow.utils;

import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Result<T> {
    T result;
    RuntimeException error;

    public Result(T result, RuntimeException error) {
        this.result = result;
        this.error = error;
    }

    public static <T> Result<T> of(Supplier<T> supplier) {
        try {
            return Result.success(supplier.get());
        } catch (RuntimeException error) {
            return Result.failure(error);
        }
    }

    public static <T> Result<T> success(T result) {
        return new Result<>(result, null);
    }

    public static <T> Result<T> failure(RuntimeException error) {
        return new Result<>(null, error);
    }

    public static <T> void throwIfFailure(Stream<Result<T>> results) {
        List<Result<T>> resultList = results.toList();
        for (Result<T> result : resultList) {
            result.getOrThrow();
        }
    }

    public T getOrThrow() {
        if (this.error != null) {
            throw this.error;
        }
        return this.result;
    }
}
